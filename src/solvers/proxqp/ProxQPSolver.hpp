#ifndef WBC_SOLVERS_PROXQP_SOLVER_HPP
#define WBC_SOLVERS_PROXQP_SOLVER_HPP

#include "../acados_qp/AcadosDenseQPSolver.hpp"

#include <string>

namespace wbc {

/**
 * @brief The ProxQPSolver class calls the dense backend of ProxQP (see https://github.com/Simple-Robotics/proxsuite)
 * through the dense QP interface of acados (DENSE_QP_PROXQP, acados/dense_qp/dense_qp_proxqp.c). This is the same
 * solver build that an acados MPC uses: the proxsuite version that acados pins, with the SIMD kernels if acados was
 * built with ACADOS_PROXQP_VECTORIZE. See AcadosDenseQPSolver for the problem formulation.
 *
 * The bounds are passed to ProxQP as box constraints, masked (absent) sides as ProxQP's infinity.
 *
 * Options, as acados maps them:
 *  - setTolerance(tol): eps_abs = tol (fields tol_stat, tol_eq, tol_ineq) and eps_duality_gap_abs = tol (fields
 *    tol_comp, tol_dual_gap). max_iter = 1000.
 *  - setMaxIter(n): max_iter = n.
 *  - setWarmStart(level): 0 is a cold start. 1 starts from the previous solution, but resets the proximal parameters
 *    (COLD_START_WITH_PREVIOUS_RESULT). 2 (default) and 3 also keep the proximal parameters
 *    (WARM_START_WITH_PREVIOUS_RESULT). acados never warm starts from a failed solve.
 *  - setDenseBackend(backend): factorization of the dense backend.
 *
 * Reference:
 * Antoine Bambade, Sarah El-Kazdadi, Adrien Taylor, Justin Carpentier. PROX-QP: Yet another Quadratic Programming Solver for
 * Robotics and beyond. RSS 2022 - Robotics: Science and Systems, Jun 2022, New York, United States. ⟨hal-03683733⟩
 */
class ProxQPSolver : public AcadosDenseQPSolver{
private:
    static QPSolverRegistry<ProxQPSolver> reg;

public:
    ProxQPSolver();

    /** Select the factorization of the dense backend: AUTOMATIC (default), PRIMAL_DUAL_LDLT or PRIMAL_LDLT (acados
      * option proxqp_dense_backend). Takes effect at the next call to solve(). Throws std::invalid_argument for an
      * unknown backend. */
    void setDenseBackend(const std::string& backend);

    /** Factorization of the dense backend, see setDenseBackend()*/
    const std::string& getDenseBackend() const { return dense_backend; }

protected:
    virtual void setSolverOptions();

private:
    std::string dense_backend = "AUTOMATIC";
};

}

#endif
