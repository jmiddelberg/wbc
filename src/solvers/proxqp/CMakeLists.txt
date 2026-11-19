SET(TARGET_NAME wbc-solvers-proxqp)

# ProxQP comes with acados (ACADOS_WITH_PROXQP) and is called through its dense QP interface, see
# ../acados_dense_qp. acados is found in ../CMakeLists.txt. acados compiles proxsuite into its own
# C wrapper (libproxsuite_c, with the SIMD kernels if ACADOS_PROXQP_VECTORIZE is ON), so wbc needs
# no proxsuite installation and no SIMD flags of its own.
wbc_require_acados_solver(ACADOS_WITH_PROXQP)

set(SOURCES ProxQPSolver.cpp)
set(HEADERS ProxQPSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-acados_dense_qp)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_link_libraries(${TARGET_NAME}
                      PUBLIC  wbc-acados_dense_qp
                      PRIVATE acados)

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/proxqp)

add_subdirectory(test)
