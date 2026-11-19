#ifndef WBC_SOLVERS_ACADOS_QP_SOLVER_HPP
#define WBC_SOLVERS_ACADOS_QP_SOLVER_HPP

#include "../../core/QPSolver.hpp"

#include <Eigen/Core>
#include <memory>
#include <string>

namespace wbc {

class HierarchicalQP;
struct QuadraticProgram;

namespace acados_qp {
class Backend;
}

/**
 * @brief Common base of the solvers that wbc calls through acados. The derived classes only select the solver and
 * its options, so that wbc needs nothing but acados to get these solvers, and uses the same solver build and the
 * same option handling as an MPC that runs on acados.
 *
 * acados has two interfaces that wbc uses, one per base class:
 *  - AcadosDenseQPSolver: the dense QP interface (DAQP, ProxQP, qpOASES)
 *  - AcadosOcpQPSolver: the OCP QP solver modules with a single stage (OSQP, Clarabel), for the solvers that acados
 *    only has as OCP QP solvers
 *
 * Both get the QuadraticProgram
 *  \f[
 *        \begin{array}{ccc}
 *        min(\mathbf{x}) & \frac{1}{2} \mathbf{x}^T\mathbf{H}\mathbf{x}+\mathbf{x}^T\mathbf{g}& \\
 *             & & \\
 *        s.t. & \mathbf{Ax} = \mathbf{b}& \\
 *             & lb(\mathbf{Cx}) \leq \mathbf{Cx} \leq ub(\mathbf{Cx})& \\
 *             & lb(\mathbf{x}) \leq \mathbf{x} \leq ub(\mathbf{x})& \\
 *        \end{array}
 *  \f]
 * Every variable of a bounded QP is a box constraint (nb = nq). A side of a bound or of a general constraint that is
 * marked absent with wbc::INF is switched off with the acados constraint mask, which every acados solver replaces by
 * its own way to drop a constraint side. The dimensions thus only change with the size of the QP, and the solver
 * memory, including the warm start, is kept from one call of solve() to the next.
 *
 * The acados memory is created again when the dimensions change, after reset() and after a change of an option.
 */
class AcadosQPSolver : public QPSolver{
public:
    virtual ~AcadosQPSolver();

    // acados keeps pointers into this object (the dimensions and the options)
    AcadosQPSolver(const AcadosQPSolver&) = delete;
    AcadosQPSolver& operator=(const AcadosQPSolver&) = delete;

    /**
     * @brief solve Solve the given quadratic program
     * @param hierarchical_qp Description of the hierarchical quadratic program to solve. Only one priority level is
     *                        implemented.
     * @param solver_output solution of the quadratic program
     * @param allow_warm_start If false, this call is cold started (acados option warm_start = 0). The next call that
     *                         allows warm starting uses the level of setWarmStart() again.
     * @throws std::runtime_error if acados does not return ACADOS_SUCCESS
     */
    virtual void solve(const HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start = true);

    /** Set tol_stat, tol_eq, tol_ineq, tol_comp and tol_dual_gap to tol and iter_max to 1000, as the MPC does
      * (model_predictive_control.cpp). Each acados solver maps these fields onto its own settings, see the derived
      * classes. Takes effect at the next call to solve(). */
    virtual bool setTolerance(double tol);

    /** Set the acados option iter_max. Overrides the limit that setTolerance() sets. Takes effect at the next call to
      * solve(). Throws std::invalid_argument for n < 1. */
    void setMaxIter(int n);

    /** Set the warm start level (acados option warm_start) that solve() uses when warm starting is allowed. 0 is a cold
      * start. The meaning of the higher levels, and which ones are valid for a QP whose matrices change in every
      * cycle, depends on the solver, see the derived classes. Throws std::invalid_argument for a level the solver
      * does not support. */
    void setWarmStart(int level);

    /** Warm start level that solve() uses when warm starting is allowed*/
    int getWarmStart() const { return warm_start; }

    /** acados return value of the last call to solve() (ACADOS_SUCCESS = 0) */
    int getReturnCode() const { return return_code; }

    /** Number of solver iterations of the last call to solve(), as acados reports it*/
    int getIterations() const;

    /** Name of the solver, as used in error messages */
    const std::string& getName() const { return name; }

protected:
    /**
     * @param backend The acados interface, see AcadosDenseQPSolver and AcadosOcpQPSolver
     * @param name Name of the solver, as used in error messages
     * @param default_warm_start Warm start level that solve() uses when warm starting is allowed
     * @param max_warm_start Highest warm start level that is valid for a QP whose matrices change in every cycle
     */
    AcadosQPSolver(std::unique_ptr<acados_qp::Backend> backend, const std::string& name, int default_warm_start, int max_warm_start);

    /** Set the solver specific options with setOption(). Called each time the acados memory is created, after the
      * generic options. */
    virtual void setSolverOptions(){}

    /** Pass an option to the acados opts_set function of the solver. Only valid inside setSolverOptions(). Note that
      * acados terminates the process on a field it does not know, so only pass validated fields and values. */
    void setOption(const std::string& field, void* value);

    /** Text for an acados return value */
    static std::string returnCodeToString(int code);

private:
    std::unique_ptr<acados_qp::Backend> backend;
    struct State;
    std::unique_ptr<State> state;   /** dimensions and bounds of the current acados memory*/

    std::string name;
    int max_warm_start;
    int warm_start;
    int applied_warm_start; /** value of the acados option warm_start in the current options, -1: none*/
    double tolerance;       /** <= 0: keep the defaults of acados*/
    int iter_max;           /** <= 0: keep the defaults of acados, or the limit that setTolerance() sets*/
    int return_code;
    bool setting_options;   /** true while the options of a new acados memory are set*/

    /** Create the acados memory for the dimensions in state and set all options */
    void create();
};

}

#endif
