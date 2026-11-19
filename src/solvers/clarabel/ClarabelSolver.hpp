#ifndef WBC_SOLVERS_CLARABEL_SOLVER_HPP
#define WBC_SOLVERS_CLARABEL_SOLVER_HPP

#include "../acados_qp/AcadosOcpQPSolver.hpp"

namespace wbc {

/**
 * @brief The ClarabelSolver class calls the interior-point conic solver Clarabel (see https://clarabel.org) through
 * the OCP QP solver module of acados (acados/ocp_qp/ocp_qp_clarabel.c). This is the same solver build and the same
 * module that an acados MPC uses with PARTIAL_CONDENSING_CLARABEL. See AcadosOcpQPSolver for the problem formulation.
 *
 * acados transcribes the QP into Clarabel's conic form: the equalities go into a ZeroCone, the sides of the bounds
 * and inequalities into a NonnegativeCone. Masked (absent) sides, and sides whose bound reaches ACADOS_INFTY (1e10),
 * are not written at all, so they do not distort the central path. The sparsity pattern is detected from the data.
 * The Clarabel solver is built at the first call and when the pattern or the set of present sides grows, otherwise its
 * data is updated in place. Presolve is off, because it would forbid these updates. Clarabel uses its default linear
 * system solver, which is qdldl in the acados build.
 *
 * Options, as acados maps them:
 *  - setTolerance(tol): tol_feas = tol (fields tol_stat, tol_eq, tol_ineq), tol_gap_abs = tol_gap_rel = tol (fields
 *    tol_comp, tol_dual_gap). max_iter = 1000.
 *  - setMaxIter(n): max_iter = n.
 *  - Warm start: Clarabel has none, so only level 0 is valid.
 *
 * acados returns ACADOS_SUCCESS only for the Clarabel status "solved". "Almost solved" (Clarabel's reduced
 * tolerances, checked after an abnormal stop) is a failure, so solve() throws.
 *
 * Reference:
 * Goulart, P.J., Chen, Y. Clarabel: An interior-point solver for conic programs with quadratic
 * objectives. https://arxiv.org/abs/2405.12762
 */
class ClarabelSolver : public AcadosOcpQPSolver{
private:
    static QPSolverRegistry<ClarabelSolver> reg;

public:
    ClarabelSolver();
};

}

#endif
