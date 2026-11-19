#ifndef WBC_SOLVERS_ACADOS_DENSE_QP_SOLVER_HPP
#define WBC_SOLVERS_ACADOS_DENSE_QP_SOLVER_HPP

#include "AcadosQPSolver.hpp"

#include <string>

namespace wbc {

/**
 * @brief Base of the solvers that wbc calls through the dense QP interface of acados
 * (acados_c/dense_qp_interface.h): DAQP, ProxQP and qpOASES. See AcadosQPSolver for the QP and the options.
 */
class AcadosDenseQPSolver : public AcadosQPSolver{
protected:
    /**
     * @param acados_solver Value of the acados enum dense_qp_solver_t, e.g. DENSE_QP_DAQP
     * @param name Name of the solver, as used in error messages
     * @param default_warm_start Warm start level that solve() uses when warm starting is allowed
     * @param max_warm_start Highest warm start level that is valid for a QP whose matrices change in every cycle
     */
    AcadosDenseQPSolver(int acados_solver, const std::string& name, int default_warm_start, int max_warm_start);
};

}

#endif
