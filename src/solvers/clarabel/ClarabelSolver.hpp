#ifndef WBC_SOLVERS_CLARABEL_SOLVER_HPP
#define WBC_SOLVERS_CLARABEL_SOLVER_HPP

#include "../../core/QPSolver.hpp"

#include <Eigen/Sparse>
#include <clarabel.hpp>

namespace wbc {

class HierarchicalQP;

/**
 * @brief The ClarabelSolver class is a wrapper for the interior-point conic solver Clarabel
 * (see https://github.com/oxfordcontrol/Clarabel.cpp). It solves problems of the shape:
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
 * Internally the problem is transcribed into Clarabel's conic form
 * \f$ min\; \frac{1}{2} x^TPx + q^Tx \; s.t. \; \hat{A}x + s = \hat{b},\; s \in \mathcal{K} \f$,
 * where the equality rows go into a ZeroCone and the (two-sided) inequality and bound rows are
 * split into pairs of NonnegativeCone rows.
 *
 * A right-hand side at or beyond getInfinity() is replaced by an actual infinity, which makes
 * Clarabel's presolve drop that row from the problem entirely. See setInfinity().
 *
 * Reference:
 * Goulart, P.J., Chen, Y. Clarabel: An interior-point solver for conic programs with quadratic
 * objectives. https://arxiv.org/abs/2405.12762
 *
 * Parameters:
 *  - See clarabel::DefaultSettings for all possible options and default values.
 */
class ClarabelSolver : public QPSolver{
private:
    static QPSolverRegistry<ClarabelSolver> reg;

public:
    ClarabelSolver();
    virtual ~ClarabelSolver(){ }

    /**
     * @brief solve Solve the given quadratic program
     * @param hierarchical_qp Description of the hierarchical quadratic program to solve. Each vector entry corresponds
     *                        to a stage in the hierarchy where the first entry has the highest priority. Currently only
     *                        one priority level is implemented.
     * @param solver_output solution of the quadratic program
     * @param allow_warm_start Unused: Clarabel is an interior-point solver and is rebuilt on every call to solve().
     */
    virtual void solve(const wbc::HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start = true);

    /** Get number of interior-point iterations of the last solve*/
    int getNter(){ return _actual_n_iter; }

    /** Set solver options. Overrides the default settings for all subsequent calls to solve().*/
    void setOptions(clarabel::DefaultSettings<double> opt){ settings = opt; }

    /**
     * @brief Set the magnitude at which a constraint bound is considered absent.
     *
     * QuadraticProgram has no way of marking a bound as unbounded, so the scenes fill unused
     * bounds with large finite sentinel values (AccelerationSceneReducedTSID for example uses
     * +/-1e4 as the default bound of every variable). For an active-set solver such a row is
     * free, but Clarabel is an interior-point method: the row stays in the KKT system that is
     * factorized in every iteration, and its distance to the iterate is orders of magnitude
     * larger than that of the rows that can actually become active, which distorts the central
     * path. Rows whose right-hand side reaches this value are therefore handed to Clarabel as
     * an actual infinity so that its presolve removes them.
     *
     * Must be larger than any bound that can genuinely become active, otherwise that constraint
     * is silently dropped. Set to infinity to disable the substitution.
     */
    void setInfinity(double value){ infinity = value; }
    double getInfinity() const { return infinity; }

protected:

    // Clarabel conic problem data (kept as members so the mutable Eigen::Ref arguments to the
    // DefaultSolver constructor bind to valid lvalues).
    Eigen::SparseMatrix<double> _P; // upper-triangular part of the (symmetric) Hessian, CSC
    Eigen::VectorXd _q;             // gradient
    Eigen::MatrixXd _A_dense;       // stacked constraint matrix (equalities + split inequalities + bounds)
    Eigen::SparseMatrix<double> _A; // sparse (CSC) view of _A_dense
    Eigen::VectorXd _b;             // stacked right-hand side

    int _actual_n_iter;

    double infinity = 1e4;          // see setInfinity()

    clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();
};

}

#endif
