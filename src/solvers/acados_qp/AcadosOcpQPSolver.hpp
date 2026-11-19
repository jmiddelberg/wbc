#ifndef WBC_SOLVERS_ACADOS_OCP_QP_SOLVER_HPP
#define WBC_SOLVERS_ACADOS_OCP_QP_SOLVER_HPP

#include "AcadosQPSolver.hpp"

#include <string>

namespace wbc {

/**
 * @brief Base of the solvers that acados only has as OCP QP solvers: OSQP and Clarabel. See AcadosQPSolver for the QP
 * and the options.
 *
 * The QuadraticProgram is an OCP QP with a single stage (N = 0) and no states: the variables are the controls of
 * stage 0, the bounds are control bounds, and the equality and inequality constraints are the general constraints
 * D u, with the equalities first. The equalities are marked as such (nge, idxe), so that a solver that has a
 * separate representation for them uses it (Clarabel: ZeroCone rows). The QP goes directly to the acados OCP QP solver
 * module (ocp_qp_osqp.c, ocp_qp_clarabel.c), without the partial condensing that an MPC puts in front of it, because
 * a single stage has nothing to condense.
 */
class AcadosOcpQPSolver : public AcadosQPSolver{
public:
    /** acados OCP QP solver modules that wbc uses*/
    enum Module{
        MODULE_OSQP,        ///< ocp_qp_osqp (acados PARTIAL_CONDENSING_OSQP without the condensing)
        MODULE_CLARABEL     ///< ocp_qp_clarabel (acados PARTIAL_CONDENSING_CLARABEL without the condensing)
    };

protected:
    /**
     * @param module The acados OCP QP solver module
     * @param name Name of the solver, as used in error messages
     * @param default_warm_start Warm start level that solve() uses when warm starting is allowed
     * @param max_warm_start Highest warm start level that is valid for a QP whose matrices change in every cycle
     */
    AcadosOcpQPSolver(Module module, const std::string& name, int default_warm_start, int max_warm_start);
};

}

#endif
