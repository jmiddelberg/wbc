#ifndef WBC_SOLVERS_ACADOS_QP_BACKEND_HPP
#define WBC_SOLVERS_ACADOS_QP_BACKEND_HPP

// Internal to the acados solvers of wbc, not installed: the interface between AcadosQPSolver and the acados
// interface that a solver uses (AcadosDenseQPSolver.cpp, AcadosOcpQPSolver.cpp).

#include <Eigen/Core>
#include <vector>

namespace wbc {

struct QuadraticProgram;

namespace acados_qp {

/** Size of a QuadraticProgram, as acados sees it */
struct Dims{
    int nv = -1;    /** variables*/
    int ne = -1;    /** equality constraints*/
    int nb = -1;    /** bounded variables: nv for a bounded QP, else 0*/
    int ng = -1;    /** inequality constraints*/

    bool operator==(const Dims& other) const { return nv == other.nv && ne == other.ne && nb == other.nb && ng == other.ng; }
    bool operator!=(const Dims& other) const { return !(*this == other); }
};

/**
 * Bounds of the variables and of the inequality constraints of a QuadraticProgram, in the layout acados expects.
 * A side that is marked absent with wbc::INF has mask 0.0. Its value is never read by the solvers, but is kept
 * finite (+-wbc::INF).
 */
struct Bounds{
    std::vector<int> idxb;  /** the bounded variables: 0..nb-1*/
    std::vector<double> lb, ub, lb_mask, ub_mask;
    std::vector<double> lg, ug, lg_mask, ug_mask;

    void resize(const Dims& dims);
    void set(const QuadraticProgram& qp);
};

/**
 * One acados interface. The life cycle is: create(), the options (setOption()), createSolver(), then any number of
 * setData() and solve(), then destroy(). destroy() is also valid on a backend that is not (completely) created.
 */
class Backend{
public:
    virtual ~Backend(){}

    /** Create the acados config, the QP data and the options for these dimensions, with the default options*/
    virtual void create(const Dims& dims) = 0;
    /** Create the solver memory. Some solvers size their memory from the options, so this comes after them.*/
    virtual void createSolver() = 0;
    /** Free everything that create() and createSolver() allocated*/
    virtual void destroy() = 0;
    /** True between create() and destroy()*/
    virtual bool created() const = 0;

    /** acados opts_set of the solver*/
    virtual void setOption(const char* field, void* value) = 0;
    /** Write the QP into the acados data, including the constraint masks*/
    virtual void setData(const QuadraticProgram& qp, const Bounds& bounds) = 0;
    /** Solve, and return the acados return value*/
    virtual int solve() = 0;
    /** Primal solution of the last solve()*/
    virtual void getSolution(Eigen::VectorXd& x) const = 0;
    /** Iterations of the last solve(), 0 if there was none*/
    virtual int iterations() const = 0;
};

}

}

#endif
