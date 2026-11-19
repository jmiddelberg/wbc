#ifndef WBC_SOLVERS_QP_OASES_SOLVER_HPP
#define WBC_SOLVERS_QP_OASES_SOLVER_HPP

#include "../acados_qp/AcadosDenseQPSolver.hpp"

namespace wbc {

/**
 * @brief The QPOASESSolver class calls qpOASES (see https://github.com/coin-or/qpOASES) through the dense QP interface
 * of acados (DENSE_QP_QPOASES, acados/dense_qp/dense_qp_qpoases.c). acados comes with qpOASES_e, the C version of
 * qpOASES, which is the same solver build that an acados MPC uses with FULL_CONDENSING_QPOASES. See
 * AcadosDenseQPSolver for the problem formulation.
 *
 * acados passes the equality constraints to qpOASES as general constraints with equal bounds (b <= Ax <= b, this
 * needs acados with commit 14b6c7f), and masked (absent) bounds as qpOASES' infinity. acados builds qpOASES_e with
 * __SUPPRESSANYOUTPUT__, so qpOASES prints nothing.
 *
 * In contrast to the earlier wrapper of the C++ version (SQProblem::hotstart), acados initializes a new QProblem in
 * every call. With warm starting, it guesses the initial working set from the dual solution of the previous call.
 *
 * Options, as acados maps them:
 *  - qpOASES runs with the MPC preset (Options_setToMPC) and terminationTolerance = 1e-4 by default.
 *  - setTolerance(tol): terminationTolerance = min(1e-4, tol). A tolerance above 1e-4 thus has no effect.
 *  - setMaxIter(n): maximum number of working set recalculations nWSR = n (default 1000).
 *  - setWarmStart(level): 0 is a cold start (QProblem_init), 1 (default) guesses the working set from the previous
 *    dual solution (QProblem_initW), without the constraints that are absent now. acados also has level 2, a hot
 *    start that is only valid for constant matrices. It is not accepted here.
 *
 * Reference:
 * Ferreau, H.J., Kirches, C., Potschka, A. et al. qpOASES: a parametric active-set algorithm for quadratic programming.
 * Math. Prog. Comp. 6, 327–363 (2014). https://doi.org/10.1007/s12532-014-0071-1
 */
class QPOASESSolver : public AcadosDenseQPSolver{
private:
    static QPSolverRegistry<QPOASESSolver> reg;

public:
    QPOASESSolver();
};

}

#endif
