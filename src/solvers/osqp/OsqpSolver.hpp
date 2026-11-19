#ifndef WBC_OSQP_SOLVER_HPP
#define WBC_OSQP_SOLVER_HPP

#include "../acados_qp/AcadosOcpQPSolver.hpp"

namespace wbc {

/**
 * @brief The OsqpSolver class calls the operator splitting solver OSQP (see https://osqp.org) through the OCP QP
 * solver module of acados (acados/ocp_qp/ocp_qp_osqp.c). This is the same solver build and the same module that an
 * acados MPC uses with PARTIAL_CONDENSING_OSQP. See AcadosOcpQPSolver for the problem formulation.
 *
 * acados writes all constraints as rows l <= A x <= u of OSQP: the equalities with l = u, the masked (absent) sides
 * as OSQP's infinity, so that OSQP treats the row as one-sided. The sparsity pattern of P and A is detected from the
 * data. The OSQP solver is set up at the first call and when the pattern grows, otherwise its data is updated in
 * place.
 *
 * Options, as acados maps them:
 *  - setTolerance(tol): eps_abs = eps_rel = tol (fields tol_stat, tol_eq, tol_ineq, tol_dual_gap). max_iter = 1000.
 *  - setMaxIter(n): max_iter = n.
 *  - setWarmStart(level): 0 starts from zero, 1 (default) starts from the primal and dual solution of the previous
 *    call (OSQP warm_starting).
 * acados changes two OSQP defaults: solution polishing is on, and termination is checked every 5 iterations.
 *
 * acados returns ACADOS_SUCCESS only for the OSQP status "solved". "Solved inaccurate" and the iteration limit are
 * failures, so solve() throws.
 *
 * Reference:
 * Stellato, B., Banjac, G., Goulart, P. et al. OSQP: an operator splitting solver for quadratic programs.
 * Math. Prog. Comp. 12, 637–672 (2020). https://doi.org/10.1007/s12532-020-00179-2
 */
class OsqpSolver : public AcadosOcpQPSolver{
private:
    static QPSolverRegistry<OsqpSolver> reg;

public:
    OsqpSolver();
};

}

#endif
