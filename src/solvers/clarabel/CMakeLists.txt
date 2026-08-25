SET(TARGET_NAME wbc-solvers-clarabel)

# Clarabel.cpp (https://github.com/oxfordcontrol/Clarabel.cpp) does not provide a CMake config
# package nor a pkg-config file, so we locate its installed headers and the (Rust-built) C library
# manually. scripts/full_install.sh installs the header tree to <prefix>/include/clarabel and
# libclarabel_c.so to <prefix>/lib. A custom location can be given via -DCLARABEL_ROOT=<prefix>.
find_package(Eigen3 3.4 CONFIG REQUIRED)

# Locations can be overridden explicitly with -DCLARABEL_INCLUDE_DIR=... -DCLARABEL_LIBRARY=...
# or pointed at an install prefix with -DCLARABEL_ROOT=<prefix>.
find_path(CLARABEL_INCLUDE_DIR
          NAMES clarabel.hpp
          HINTS ${CLARABEL_ROOT} ENV CLARABEL_ROOT
          PATHS /usr/local /usr /opt/clarabel
          PATH_SUFFIXES clarabel include/clarabel)
find_library(CLARABEL_LIBRARY
             NAMES clarabel_c libclarabel_c
             HINTS ${CLARABEL_ROOT} ENV CLARABEL_ROOT
             PATHS /usr/local /usr /opt/clarabel
             PATH_SUFFIXES lib lib64)

if(NOT CLARABEL_INCLUDE_DIR OR NOT CLARABEL_LIBRARY)
    message(FATAL_ERROR
        "Clarabel not found:\n"
        "    CLARABEL_INCLUDE_DIR = ${CLARABEL_INCLUDE_DIR}   (searched for clarabel.hpp)\n"
        "    CLARABEL_LIBRARY     = ${CLARABEL_LIBRARY}   (searched for libclarabel_c)\n"
        "Build Clarabel.cpp and install its 'include' tree to <prefix>/include/clarabel and\n"
        "libclarabel_c.so to <prefix>/lib BEFORE configuring wbc (see docker/full_install steps),\n"
        "or pass -DCLARABEL_ROOT=<prefix>  (or -DCLARABEL_INCLUDE_DIR=... -DCLARABEL_LIBRARY=...).")
endif()
message(STATUS "Found Clarabel: ${CLARABEL_LIBRARY} (headers: ${CLARABEL_INCLUDE_DIR})")

# Which linear system solvers Clarabel offers, and the layout of its settings struct, depend on the
# cargo features libclarabel_c was built with, so the C++ side has to be compiled with the matching
# FEATURE_* macros (Clarabel.cpp's own CMakeLists does this for in-tree builds). Since we consume a
# pre-built library, the features have to be declared here. The defaults match a plain
# `cargo build --release` of the rust wrapper, which provides the qdldl solver only.
option(CLARABEL_FEATURE_FAER_SPARSE   "libclarabel_c was built with the faer-sparse cargo feature (adds the 'faer' linear system solver)"           OFF)
option(CLARABEL_FEATURE_PARDISO_MKL   "libclarabel_c was built with the pardiso-mkl cargo feature (adds the 'pardiso-mkl' linear system solver)"     OFF)
option(CLARABEL_FEATURE_PARDISO_PANUA "libclarabel_c was built with the pardiso-panua cargo feature (adds the 'pardiso-panua' linear system solver)" OFF)
option(CLARABEL_FEATURE_SDP           "libclarabel_c was built with one of the sdp cargo features"                                                   OFF)

# A mismatch between the options above and the library silently changes the meaning of clarabel's
# settings struct, so warn about the cases that can be told apart by looking at the library.
function(check_clarabel_feature option pattern feature)
    file(STRINGS ${CLARABEL_LIBRARY} matches REGEX "${pattern}" LIMIT_COUNT 1)
    if(${option} AND NOT matches)
        message(WARNING
            "-D${option}=ON was given, but ${CLARABEL_LIBRARY} contains no trace of the '${feature}' "
            "cargo feature. If it really was not built with that feature, the settings passed to Clarabel "
            "will be misinterpreted. Rebuild it with 'cargo build --release --features clarabel/${feature}'.")
    elseif(NOT ${option} AND matches)
        message(WARNING
            "${CLARABEL_LIBRARY} looks like it was built with the '${feature}' cargo feature, but "
            "${option} is OFF, so the settings passed to Clarabel may be misinterpreted and the linear "
            "system solvers of that feature are not selectable. Pass -D${option}=ON.")
    endif()
endfunction()

# The Pardiso variants share the pardiso-wrapper crate, so they cannot be told apart by their symbols
check_clarabel_feature(CLARABEL_FEATURE_FAER_SPARSE "faer_ldl" "faer-sparse")
if(CLARABEL_FEATURE_PARDISO_MKL OR CLARABEL_FEATURE_PARDISO_PANUA)
    check_clarabel_feature(CLARABEL_FEATURE_PARDISO_MKL "pardiso" "pardiso-mkl")
else()
    check_clarabel_feature(CLARABEL_FEATURE_PARDISO_PANUA "pardiso" "pardiso-panua")
endif()
check_clarabel_feature(CLARABEL_FEATURE_SDP "chordal" "sdp")

set(CLARABEL_FEATURE_DEFINITIONS "")
if(CLARABEL_FEATURE_FAER_SPARSE)
    list(APPEND CLARABEL_FEATURE_DEFINITIONS FEATURE_FAER_SPARSE)
endif()
if(CLARABEL_FEATURE_PARDISO_MKL)
    list(APPEND CLARABEL_FEATURE_DEFINITIONS FEATURE_PARDISO_MKL)
endif()
if(CLARABEL_FEATURE_PARDISO_PANUA)
    list(APPEND CLARABEL_FEATURE_DEFINITIONS FEATURE_PARDISO_PANUA)
endif()
if(CLARABEL_FEATURE_PARDISO_MKL OR CLARABEL_FEATURE_PARDISO_PANUA)
    list(APPEND CLARABEL_FEATURE_DEFINITIONS FEATURE_PARDISO_ANY)
endif()
if(CLARABEL_FEATURE_SDP)
    list(APPEND CLARABEL_FEATURE_DEFINITIONS FEATURE_SDP)
endif()
if(CLARABEL_FEATURE_DEFINITIONS)
    message(STATUS "Clarabel features: ${CLARABEL_FEATURE_DEFINITIONS}")
endif()

set(SOURCES ClarabelSolver.cpp)
set(HEADERS ClarabelSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

# Anything that includes ClarabelSolver.hpp pulls in clarabel.hpp and therefore needs the same
# feature macros as this library
foreach(definition ${CLARABEL_FEATURE_DEFINITIONS})
    list(APPEND PKGCONFIG_CFLAGS -D${definition})
endforeach()
string (REPLACE ";" " " PKGCONFIG_CFLAGS "${PKGCONFIG_CFLAGS}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_include_directories(${TARGET_NAME} PUBLIC ${CLARABEL_INCLUDE_DIR})
target_compile_definitions(${TARGET_NAME} PUBLIC ${CLARABEL_FEATURE_DEFINITIONS})
target_link_libraries(${TARGET_NAME} PUBLIC
                      wbc-core
                      Eigen3::Eigen
                      ${CLARABEL_LIBRARY})
# libclarabel_c is a Rust cdylib and needs these system libraries at link time.
target_link_libraries(${TARGET_NAME} PUBLIC pthread dl m)

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

# The Clarabel C++ interface requires C++17; the rest of wbc builds with C++14.
target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/clarabel)

add_subdirectory(test)
