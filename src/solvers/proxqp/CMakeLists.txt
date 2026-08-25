SET(TARGET_NAME wbc-solvers-proxqp)

# proxsuite ships a single proxsuite.pc, whether or not it was built with
# BUILD_WITH_VECTORIZATION_SUPPORT, so its vectorized variant is only reachable
# through the CMake package - hence find_package rather than pkg_search_module.
find_package(proxsuite REQUIRED)

set(SOURCES ProxQPSolver.cpp)
set(HEADERS ProxQPSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES proxsuite)
list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

# ProxQP's hand-vectorized dense kernels, worth roughly 25-30% on QPs of the size
# the scenes produce. proxsuite::proxsuite-vectorized carries only the
# PROXSUITE_VECTORIZE define and simde; simde emits real SIMD instructions only if
# the compiler additionally targets an ISA that has them, so without the
# architecture flags below the define gains nothing at all.
#
# Both are applied PRIVATEly so they do not reach wbc's consumers. proxsuite is
# header-only, which means every translation unit that instantiates its templates
# has to agree on these flags: instantiations that differ are an ODR violation and
# the dynamic linker silently picks one of them for all callers. ProxQPSolver.hpp
# does include a proxsuite header, so any new translation unit that includes it
# must link ${WBC_PROXQP_VECTORIZED} and add ${WBC_PROXQP_ARCH_FLAGS} too - see
# test/CMakeLists.txt for the one place that already does.
set(WBC_PROXQP_VECTORIZED "")
set(WBC_PROXQP_ARCH_FLAGS "")
if(TARGET proxsuite::proxsuite-vectorized)
    set(WBC_PROXQP_VECTORIZED proxsuite::proxsuite-vectorized)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        set(WBC_PROXQP_ARCH_FLAGS -mavx2 -mfma)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        set(WBC_PROXQP_ARCH_FLAGS "")   # NEON is part of the ARMv8-A baseline
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        set(WBC_PROXQP_ARCH_FLAGS -mfpu=neon)
    endif()
    message(STATUS "wbc-solvers-proxqp: ProxQP SIMD kernels enabled (${WBC_PROXQP_ARCH_FLAGS})")
else()
    message(STATUS "wbc-solvers-proxqp: proxsuite built without BUILD_WITH_VECTORIZATION_SUPPORT, using the scalar kernels")
endif()

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_link_libraries(${TARGET_NAME}
                      PUBLIC  wbc-core
                              proxsuite::proxsuite
                      PRIVATE ${WBC_PROXQP_VECTORIZED})
target_compile_options(${TARGET_NAME} PRIVATE ${WBC_PROXQP_ARCH_FLAGS})

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/proxqp)

add_subdirectory(test)
