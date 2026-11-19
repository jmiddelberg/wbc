#ifndef WBC_SOLVERS_DAQP_SOLVER_HPP
#define WBC_SOLVERS_DAQP_SOLVER_HPP

#include "../acados_dense_qp/AcadosDenseQPSolver.hpp"

namespace wbc {

/**
 * @brief The DAQPSolver class calls the dual active-set solver DAQP (see https://github.com/darnstrom/daqp)
 * through the dense QP interface of acados (DENSE_QP_DAQP, acados/dense_qp/dense_qp_daqp.c). This is the same
 * solver build and the same interface that an acados MPC uses with FULL_CONDENSING_DAQP. See AcadosDenseQPSolver
 * for the problem formulation.
 *
 * acados transforms the QP into the least distance problem that DAQP solves. Equality constraints stay in the
 * working set at all times, and masked (absent) bounds are set to DAQP's infinity, so they never enter it.
 *
 * Options, as acados maps them:
 *  - setTolerance(tol): primal_tol = dual_tol = tol (field tol_ineq). The other tolerances are enforced by the
 *    working set itself. iter_limit = 1000.
 *  - setMaxIter(n): iter_limit = n.
 *  - setWarmStart(level): 0 starts from an empty working set, 1 (default) starts from the working set of the previous
 *    solution, without the constraints that are absent now. acados also has level 2, which keeps the factorization
 *    of the previous call and is thus only valid for constant matrices. It is not accepted here.
 *
 * Reference:
 * Arnström, D., Bemporad, A., Axehill, D. A Dual Active-Set Solver for Embedded Quadratic
 * Programming Using Recursive LDL' Updates. IEEE Transactions on Automatic Control, 2022.
 * https://doi.org/10.1109/TAC.2022.3176430
 */
class DAQPSolver : public AcadosDenseQPSolver{
private:
    static QPSolverRegistry<DAQPSolver> reg;

public:
    DAQPSolver();
};

}

#endif // WBC_SOLVERS_DAQP_SOLVER_HPP
