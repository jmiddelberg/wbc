#ifndef WBC_SOLVERS_CLARABEL_SOLVER_HPP
#define WBC_SOLVERS_CLARABEL_SOLVER_HPP

#include "../../core/QPSolver.hpp"
#include "../../core/QuadraticProgram.hpp"

#include <Eigen/Sparse>
#include <clarabel.hpp>
#include <string>
#include <vector>

namespace wbc {

class HierarchicalQP;

/**
 * @brief Linear system solver that Clarabel uses to factorize the KKT system in every
 * interior-point iteration, see https://clarabel.org/stable/user_guide_linsolvers/.
 *
 * Only qdldl (and auto, which resolves to qdldl unless faer is compiled in) is always available.
 * The other methods have to be compiled into libclarabel_c via a cargo feature, and wbc has to be
 * configured with the matching cmake option, because both the numbering of Clarabel's solver enum
 * and the layout of its settings struct depend on which features are enabled:
 *
 *   | linear solver                    | cargo feature | wbc cmake option                    |
 *   | clarabel_linsolver_faer          | faer-sparse   | -DCLARABEL_FEATURE_FAER_SPARSE=ON   |
 *   | clarabel_linsolver_pardiso_mkl   | pardiso-mkl   | -DCLARABEL_FEATURE_PARDISO_MKL=ON   |
 *   | clarabel_linsolver_pardiso_panua | pardiso-panua | -DCLARABEL_FEATURE_PARDISO_PANUA=ON |
 *
 * The Pardiso methods additionally require the corresponding Pardiso library at runtime
 * (libmkl_rt.so via MKLROOT/MKL_PARDISO_PATH, libpardiso.so via PARDISO_PATH), and Panua Pardiso
 * requires a commercial license. See ClarabelSolver::availableLinearSolvers() to query what the
 * current build supports.
 */
enum ClarabelLinearSolver{
    clarabel_linsolver_auto = 0,      ///< Let Clarabel choose the method (the Clarabel default)
    clarabel_linsolver_qdldl,         ///< QDLDL, the LDL' factorization built into Clarabel
    clarabel_linsolver_faer,          ///< Faer supernodal LDL' factorization, multi-threaded
    clarabel_linsolver_pardiso_mkl,   ///< Pardiso from the Intel MKL (x86_64 only)
    clarabel_linsolver_pardiso_panua  ///< Panua Pardiso (commercial license required)
};

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
 *  - The linear system solver used for the KKT system can be selected with setLinearSolver(),
 *    see ClarabelLinearSolver.
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

    /**
     * @brief Set solver options. Overrides the default settings for all subsequent calls to
     * solve(). Also overrides the linear system solver selected by setLinearSolver() with
     * opt.direct_solve_method.
     */
    void setOptions(clarabel::DefaultSettings<double> opt);

    /** Return the current solver options*/
    clarabel::DefaultSettings<double> getOptions() const { return settings; }

    /**
     * @brief Select the linear system solver that Clarabel uses to factorize the KKT system in
     * every interior-point iteration, see ClarabelLinearSolver.
     *
     * The default is clarabel_linsolver_qdldl rather than Clarabel's own default
     * clarabel_linsolver_auto: the QPs that the wbc scenes produce have dense H and C, for which
     * qdldl was measured to be the fastest of the available methods, whereas auto switches to the
     * supernodal faer factorization from roughly 50 variables on.
     *
     * @param method The method to use. Throws std::runtime_error if it is not available in this
     *               build, see linearSolverAvailable().
     */
    void setLinearSolver(ClarabelLinearSolver method);

    /**
     * @brief Select the linear system solver by name, using the names of the Clarabel
     * documentation: "auto", "qdldl", "faer", "pardiso-mkl" or "pardiso-panua". Throws
     * std::invalid_argument if the name is unknown and std::runtime_error if the method is not
     * available in this build.
     */
    void setLinearSolver(const std::string& method);

    /** Get the linear system solver that subsequent calls to solve() will request*/
    ClarabelLinearSolver getLinearSolver() const { return linear_solver; }

    /**
     * @brief Get the linear system solver that the last call to solve() actually used, as
     * reported by Clarabel. This is the method that clarabel_linsolver_auto resolved to, and it
     * is only meaningful after solve() has been called at least once.
     */
    ClarabelLinearSolver getLinearSolverUsed() const { return _linsolver_used; }

    /** Check whether the given linear system solver is compiled into this build, see ClarabelLinearSolver*/
    static bool linearSolverAvailable(ClarabelLinearSolver method);

    /** Names of all linear system solvers that are compiled into this build*/
    static std::vector<std::string> availableLinearSolvers();

    /** Name of the given linear system solver, as used in the Clarabel documentation*/
    static std::string linearSolverName(ClarabelLinearSolver method);

    /** Inverse of linearSolverName(). Throws std::invalid_argument if the name is unknown*/
    static ClarabelLinearSolver linearSolverFromName(const std::string& name);

    /**
     * @brief Set the magnitude at which a constraint bound is considered absent.
     *
     * Defaults to wbc::INF, the value the scenes use to mark a bound as absent. For an
     * active-set solver such a row is free, but Clarabel is an interior-point method: the row
     * stays in the KKT system that is factorized in every iteration, and its distance to the
     * iterate is orders of magnitude larger than that of the rows that can actually become
     * active, which distorts the central path. Rows whose right-hand side reaches this value are
     * therefore handed to Clarabel as an actual infinity so that its presolve removes them.
     *
     * Lowering it below wbc::INF also discards bounds that are merely large but were meant to be
     * enforced, so only do that if you know that no bound of that magnitude can become active.
     * Set it to infinity to disable the substitution altogether.
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
    ClarabelLinearSolver _linsolver_used = clarabel_linsolver_auto; // method the last solve() used

    double infinity = INF;          // see setInfinity()

    ClarabelLinearSolver linear_solver = clarabel_linsolver_qdldl;  // see setLinearSolver()

    clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();
};

}

#endif
