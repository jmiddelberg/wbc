#ifndef WBC_SOLVERS_DAQP_SOLVER_HPP
#define WBC_SOLVERS_DAQP_SOLVER_HPP

#include "../../core/QPSolver.hpp"

#include <Eigen/Core>
#include <daqp/api.h>

namespace wbc {

class HierarchicalQP;

/**
 * @brief The DAQPSolver class is a wrapper for the dual active-set solver DAQP
 * (see https://github.com/darnstrom/daqp). It solves problems of the shape:
 *  \f[
 *        \begin{array}{ccc}
 *        min(\mathbf{x}) & \frac{1}{2} \mathbf{x}^T\mathbf{H}\mathbf{x}+\mathbf{x}^T\mathbf{g}& \\
 *             & & \\
 *        s.t. & \mathbf{Ax} = \mathbf{b}& \\
 *             & lb(\mathbf{Cx}) \leq \mathbf{Cx} \leq ub(\mathbf{Cx})& \\
 *             & lb(\mathbf{x}) \leq \mathbf{x} \leq ub(\mathbf{x})& \\
 *        \end{array}
 *  \f]
 *
 * DAQP expects all constraints in the single two-sided form
 * \f$ b_l \leq [I;\hat{A}]x \leq b_u \f$, where the first ms rows are the simple bounds on the
 * first ms variables and the remaining rows are the (row-major) general constraint matrix.
 * Equality rows are passed with identical bounds and marked as active and immutable
 * (sense = DAQP_ACTIVE|DAQP_IMMUTABLE), which is how DAQP encodes equality constraints.
 *
 * Since DAQP is an active-set method, subsequent calls to solve() re-use the solver workspace
 * and, if warm starting is allowed, the active set of the previous solution. The workspace is
 * rebuilt whenever the problem dimensions change or warm starting is switched off.
 *
 * Reference:
 * Arnström, D., Bemporad, A., Axehill, D. A Dual Active-Set Solver for Embedded Quadratic
 * Programming Using Recursive LDL' Updates. IEEE Transactions on Automatic Control, 2022.
 * https://doi.org/10.1109/TAC.2022.3176430
 *
 * Parameters:
 *  - See DAQPSettings in daqp/types.h for all possible options,
 *    https://darnstrom.github.io/daqp/parameters for their meaning and default values.
 */
class DAQPSolver : public QPSolver{
private:
    static QPSolverRegistry<DAQPSolver> reg;

public:
    DAQPSolver();
    virtual ~DAQPSolver();

    // The solver owns a DAQP workspace, which holds pointers into the members of this class
    DAQPSolver(const DAQPSolver&) = delete;
    DAQPSolver& operator=(const DAQPSolver&) = delete;

    /**
     * @brief solve Solve the given quadratic program
     * @param hierarchical_qp Description of the hierarchical quadratic program to solve. Each vector entry corresponds
     *                        to a stage in the hierarchy where the first entry has the highest priority. Currently only
     *                        one priority level is implemented.
     * @param solver_output solution of the quadratic program
     * @param allow_warm_start If false, the solver workspace is set up from scratch and the active set of the previous
     *                         solution is discarded.
     */
    virtual void solve(const wbc::HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start = true);

    /** Get number of active set iterations of the last solve*/
    int getNter(){ return _result.iter; }

    /** Get the exit flag of the last solve (see the DAQP_EXIT_* defines in daqp/constants.h)*/
    int getExitFlag(){ return _result.exitflag; }

    /** Return the current solver options*/
    DAQPSettings getOptions(){ return _settings; }

    /** Set solver options. Enforces a reconfiguration at the next call to solve().*/
    void setOptions(const DAQPSettings& opt){ _settings = opt; reset(); }

protected:

    /** Free the DAQP workspace (if any). Does not touch the solver options.*/
    void freeWorkspace();

    DAQPWorkspace _work;
    DAQPSettings _settings;
    DAQPResult _result;
    DAQPProblem _problem;

    // Problem data in the layout expected by DAQP: the constraint matrix stacks the equalities on
    // top of the inequalities and is stored row-major, the bound vectors are the simple bounds
    // followed by the general constraints:
    //     _A = [A; C],  _b_upper = [upper_x; b; upper_y],  _b_lower = [lower_x; b; lower_y]
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _A;
    Eigen::VectorXd _b_upper;
    Eigen::VectorXd _b_lower;
    Eigen::VectorXi _sense;   // Constraint types (equality/inequality, initial active set)

    Eigen::VectorXd _x;       // Primal solution, written by DAQP
    Eigen::VectorXd _lam;     // Dual solution, written by DAQP

    // Dimensions the workspace is currently set up for
    int _n_var;               // Number of variables
    int _n_cstr;              // Total number of constraints (simple bounds + general)
    int _n_bnd;               // Number of simple bounds
};

}

#endif // WBC_SOLVERS_DAQP_SOLVER_HPP
