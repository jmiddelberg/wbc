SET(TARGET_NAME wbc-solvers-daqp)

# DAQP (https://github.com/darnstrom/daqp) installs a CMake config package to <prefix>/cmake, which
# is where find_package looks first. If that package is not found, e.g. because only the headers and
# the library were copied into place, fall back to locating those manually. Both lookups can be
# pointed at a custom install prefix with -DDAQP_ROOT=<prefix>.
find_package(daqp CONFIG QUIET HINTS ${DAQP_ROOT} ENV DAQP_ROOT)

if(daqp_FOUND)
    message(STATUS "Found DAQP: ${daqp_DIR} (CMake package)")
    set(DAQP_LIBRARIES daqp)
    set(DAQP_INCLUDE_DIRS "")
else()
    find_path(DAQP_INCLUDE_DIR
              NAMES daqp/api.h
              HINTS ${DAQP_ROOT} ENV DAQP_ROOT
              PATHS /usr/local /usr
              PATH_SUFFIXES include)
    find_library(DAQP_LIBRARY
                 NAMES daqp
                 HINTS ${DAQP_ROOT} ENV DAQP_ROOT
                 PATHS /usr/local /usr
                 PATH_SUFFIXES lib lib64)

    if(NOT DAQP_INCLUDE_DIR OR NOT DAQP_LIBRARY)
        message(FATAL_ERROR
            "DAQP not found:\n"
            "    DAQP_INCLUDE_DIR = ${DAQP_INCLUDE_DIR}   (searched for daqp/api.h)\n"
            "    DAQP_LIBRARY     = ${DAQP_LIBRARY}   (searched for libdaqp)\n"
            "Build and install DAQP BEFORE configuring wbc (see scripts/full_install.sh), or pass\n"
            "-DDAQP_ROOT=<prefix>  (or -DDAQP_INCLUDE_DIR=... -DDAQP_LIBRARY=...).")
    endif()

    message(STATUS "Found DAQP: ${DAQP_LIBRARY} (headers: ${DAQP_INCLUDE_DIR})")
    set(DAQP_LIBRARIES ${DAQP_LIBRARY})
    set(DAQP_INCLUDE_DIRS ${DAQP_INCLUDE_DIR})
endif()

set(SOURCES DAQPSolver.cpp)
set(HEADERS DAQPSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_include_directories(${TARGET_NAME} PUBLIC ${DAQP_INCLUDE_DIRS})
target_link_libraries(${TARGET_NAME} PUBLIC
                      wbc-core
                      ${DAQP_LIBRARIES}
                      m)

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/daqp)

add_subdirectory(test)
