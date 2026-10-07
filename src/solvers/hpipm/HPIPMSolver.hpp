#ifndef WBC_SOLVERS_HPIPM_SOLVER_HPP
#define WBC_SOLVERS_HPIPM_SOLVER_HPP

#include "../../core/QPSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include <acados_c/dense_qp_interface.h>
#include <string>

namespace wbc {

class HPIPMSolver : public QPSolver{
private:
    static QPSolverRegistry<HPIPMSolver> reg;
    std::vector<int> idxb;             /** Indices of the variables that actually carry a bound */
    std::vector<double> lb, ub;        /** Bounds of those variables, in idxb order */
    std::vector<double> lb_mask, ub_mask, lg_mask, ug_mask; /** 1.0 = side is enforced, 0.0 = absent */
    dense_qp_in *qp_in;
    dense_qp_dims dims;
    dense_qp_out *qp_out;
    void *opts;
    dense_qp_solver_plan plan;
    dense_qp_solver *qp_solver;
    qp_solver_config *config;

    /** Collect the variables that carry a bound on at least one side into idxb, and fill the
      * bound and mask vectors. Returns true if the set of bounded variables changed, in which
      * case the HPIPM problem has to be created again for the new dimensions. */
    bool updateBounds(const wbc::QuadraticProgram &qp);

    std::string returnCodeToString(int code);

public:
    HPIPMSolver();
    virtual ~HPIPMSolver();

    /**
     * @brief solve Solve the given quadratic program
     * @param hierarchical_qp Description of the hierarchical quadratic program to solve.
     * @param solver_output solution of the quadratic program
     */
    virtual void solve(const wbc::HierarchicalQP &hierarchical_qp, Eigen::VectorXd &solver_output, bool allow_warm_start = true);

    /** Set tol_stat, tol_eq, tol_ineq, tol_comp and tol_dual_gap to tol, as the MPC does, and raise iter_max to 1000*/
    virtual bool setTolerance(double tol);

    /** Select the HPIPM mode: SPEED_ABS, SPEED (default), BALANCE or ROBUST*/
    virtual bool setMode(const std::string& mode);

private:
    std::string mode = "SPEED";
    double tolerance = -1.0;           /** <= 0: keep the tolerances of the mode preset */
};
}

#endif
