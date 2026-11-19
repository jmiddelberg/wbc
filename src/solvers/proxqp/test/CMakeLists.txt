add_executable(test_proxqp_solver test_proxqp_solver.cpp)
target_link_libraries(test_proxqp_solver
                      wbc-solvers-proxqp
                      Boost::unit_test_framework)

add_test(NAME test_proxqp_solver COMMAND test_proxqp_solver)
