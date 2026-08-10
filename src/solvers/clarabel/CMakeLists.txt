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

set(SOURCES ClarabelSolver.cpp)
set(HEADERS ClarabelSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_include_directories(${TARGET_NAME} PUBLIC ${CLARABEL_INCLUDE_DIR})
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
