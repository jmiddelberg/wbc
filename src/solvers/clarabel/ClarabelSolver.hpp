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

protected:

    // Clarabel conic problem data (kept as members so the mutable Eigen::Ref arguments to the
    // DefaultSolver constructor bind to valid lvalues).
    Eigen::SparseMatrix<double> _P; // upper-triangular part of the (symmetric) Hessian, CSC
    Eigen::VectorXd _q;             // gradient
    Eigen::MatrixXd _A_dense;       // stacked constraint matrix (equalities + split inequalities + bounds)
    Eigen::SparseMatrix<double> _A; // sparse (CSC) view of _A_dense
    Eigen::VectorXd _b;             // stacked right-hand side

    int _actual_n_iter;

    clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();
};

}

#endif
