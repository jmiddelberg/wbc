#ifndef WBC_CORE_QUADRATIC_PROGRAM_HPP
#define WBC_CORE_QUADRATIC_PROGRAM_HPP

#include <Eigen/Core>
#include <vector>

namespace wbc{

/**
 * @brief Magnitude that marks a constraint bound as absent (unbounded).
 *
 * QuadraticProgram stores every bound as a plain double and has no separate way of saying that a
 * constraint is one-sided or that a variable is free. Scenes and constraints therefore fill unused
 * bounds with this value instead of an arbitrary large number, so that the solver wrappers can
 * recognise them and hand them to their solver in whatever way that solver expects: some accept an
 * infinity of their own (OSQP, qpOASES, ProxQP, DAQP), some remove the row (Clarabel, qpSWIFT,
 * eiquadprog) and some switch the constraint off with a mask (HPIPM).
 *
 * A large finite value is used rather than std::numeric_limits<double>::infinity() on purpose: code
 * that does not know about the convention degrades to a very loose but finite bound, whereas an
 * actual infinity would silently turn into NaN inside most solvers.
 *
 * Use isUnbounded() / hasLowerBound() / hasUpperBound() to test for it, never an equality
 * comparison - a bound is considered absent as soon as it reaches this magnitude.
 */
constexpr double INF = 1e20;

/** True if the given lower bound is absent, i.e. the constraint is not bounded from below */
inline bool hasLowerBound(double lower){ return lower > -INF; }

/** True if the given upper bound is absent, i.e. the constraint is not bounded from above */
inline bool hasUpperBound(double upper){ return upper < INF; }

/** True if the given bound is absent on the side it is used for */
inline bool isUnbounded(double bound){ return bound <= -INF || bound >= INF; }

/**
 * @brief Describes a quadratic program of the form
 *  \f[
 *        \begin{array}{ccc}
 *        min(\mathbf{x}) & \frac{1}{2} \mathbf{x}^T\mathbf{H}\mathbf{x}+\mathbf{x}^T\mathbf{g}& \\
 *             & & \\
 *        s.t. & lb(\mathbf{Ax}) \leq \mathbf{Ax} \leq ub(\mathbf{Ax})& \\
 *             & lb(\mathbf{x}) \leq \mathbf{x} \leq ub(\mathbf{x})& \\
 *        \end{array}
 *  \f]
 */
struct QuadraticProgram{

    int nq;                 /** Number of variables */
    int neq;                /** Number of equalities constraints for this prio*/
    int nin;                /** Number of inequalities constraints for this prio*/

    bool bounded;            /** Contains simple boiunds for the variables */

    Eigen::MatrixXd H;       /** Hessian Matrix (nq x nq) */
    Eigen::VectorXd g;       /** Gradient vector (nq x 1) */
    Eigen::MatrixXd A;       /** Equalities constraint matrix (neq x nq) */
    Eigen::VectorXd b;       /** Equalities constraint vector (neq x 1) */
    Eigen::MatrixXd C;       /** Inequalities constraint matrix (nin x nq) */
    Eigen::VectorXd lower_y; /** Lower bound of the constraint vector (nin x 1) */
    Eigen::VectorXd upper_y; /** Upper bound of the constraint vector (nin x 1) */
    Eigen::VectorXd lower_x; /** Lower bound of the solution vector (nq x 1) */
    Eigen::VectorXd upper_x; /** Upper bound of the solution vector (nq x 1) */
    Eigen::VectorXd Wy;      /** Constraint weights (nc x 1). Default entry is 1. */

    /** Initialize all variables with NaN */
    void resize(uint nq, uint neq, uint nin, bool bounds);

    /** Check if matrix and vectors dims match with nq, neq, nin. Throw exception if not **/
    bool isValid() const;

    /** Print content to console*/
    void print() const;

};

/**
 * @brief Describes a hierarchy of quadratic programs
 */
struct HierarchicalQP{
    std::vector<QuadraticProgram> prios;           /** Hierarchical organized QPs. The first entriy is the highest priority*/
    Eigen::VectorXd Wq;                            /** Joint weights (all joints) */

    size_t size() const {
        return prios.size();
    }
    QuadraticProgram& operator[](int i) {
        return prios[i];
    }
    const QuadraticProgram& operator[](int i) const {
        return prios[i];
    }
    void operator<<(QuadraticProgram& qp) {
        prios.push_back(qp);
    }
    void resize(const size_t &n){prios.resize(n);}
};

}

#endif // WBC_CORE_QUADRATIC_PROGRAM_HPP
