add_executable(test_clarabel_solver test_clarabel_solver.cpp)
target_link_libraries(test_clarabel_solver
                      wbc-solvers-clarabel
                      Boost::unit_test_framework)

add_test(NAME test_clarabel_solver COMMAND test_clarabel_solver)
