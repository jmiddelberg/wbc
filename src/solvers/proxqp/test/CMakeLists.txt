add_executable(test_proxqp_solver test_proxqp_solver.cpp)
target_link_libraries(test_proxqp_solver
                      wbc-solvers-proxqp
                      ${WBC_PROXQP_VECTORIZED}
                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
# test_proxqp_solver.cpp includes ProxQPSolver.hpp, so it instantiates proxsuite
# templates itself and must be built with the same vectorization flags as the
# solver - otherwise the two copies differ, see ../CMakeLists.txt
target_compile_options(test_proxqp_solver PRIVATE ${WBC_PROXQP_ARCH_FLAGS})

add_test(NAME test_proxqp_solver COMMAND test_proxqp_solver)
