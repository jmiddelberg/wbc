# DAQP, ProxQP and qpOASES are called through the dense QP interface of acados, HPIPM as well. OSQP
# and Clarabel are called through the OCP QP solver modules of acados. qpOASES is always built, so
# wbc needs acados. acados installs a CMake package to <prefix>/cmake, which find_package finds for a
# prefix like /opt/acados. Otherwise pass -Dacados_DIR=<prefix>/cmake.
find_package(acados REQUIRED)

# acados terminates the process at runtime (dense_qp_config_create) or does not link (OCP QP solver
# modules) if it was built without the solver that is requested, so this is checked when wbc is
# configured
function(wbc_require_acados_solver option)
    if(NOT ${option})
        message(FATAL_ERROR "acados (${acados_DIR}) was built without ${option}. Build acados with "
                            "-D${option}=ON, or switch the corresponding wbc solver off.")
    endif()
endfunction()

add_subdirectory(acados_qp)
add_subdirectory(qpoases)
if(SOLVER_EIQUADPROG)
    add_subdirectory(eiquadprog)
endif()
if(SOLVER_QPSWIFT)
    add_subdirectory(qpswift)
endif()
if(SOLVER_PROXQP)
    add_subdirectory(proxqp)
endif()
if(SOLVER_OSQP)
    add_subdirectory(osqp)
endif()
if(SOLVER_HPIPM)
    add_subdirectory(hpipm)
endif()
if(SOLVER_CLARABEL)
    add_subdirectory(clarabel)
endif()
if(SOLVER_DAQP)
    add_subdirectory(daqp)
endif()
