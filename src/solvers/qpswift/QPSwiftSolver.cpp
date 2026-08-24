#include "QPSwiftSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include "../../tools/Logger.hpp"

namespace wbc {

QPSolverRegistry<QPSwiftSolver> QPSwiftSolver::reg("qpswift");

QPSwiftSolver::QPSwiftSolver(){
    my_qp = 0;
    // toQpSwift() unconditionally copies these onto the freshly created problem, so they have to
    // hold qpSWIFT's own defaults rather than whatever was on the stack: an uninitialized maxit
    // and verbose make the solver iterate (and print) without end.
    options.maxit  = MAXIT;
    options.reltol = RELTOL;
    options.abstol = ABSTOL;
    options.sigma  = SIGMA;
    options.verbose = VERBOSE;
}

QPSwiftSolver::~QPSwiftSolver(){
    if(my_qp)
        QP_CLEANUP_dense(my_qp);
}

/// qpSWIFT expects one-sided inequalities Gx <= h, so every two-sided constraint of the wbc problem
/// has to be split into two rows. Sides marked with wbc::INF are simply not emitted: qpSWIFT is an
/// interior-point method, so an unbounded row is not free - it stays in the KKT system that is
/// factorized in every iteration and its barrier term pulls the iterate towards a bound that can
/// never become active. Note that this makes the number of inequality rows problem dependent, which
/// is why countInequalities() is evaluated on every call.
uint QPSwiftSolver::countInequalities(const wbc::QuadraticProgram &qp){
    uint n = 0;
    for(int i = 0; i < qp.nin; i++){
        if(hasUpperBound(qp.upper_y(i))) n++;
        if(hasLowerBound(qp.lower_y(i))) n++;
    }
    if(qp.bounded){
        for(int i = 0; i < qp.nq; i++){
            if(hasUpperBound(qp.upper_x(i))) n++;
            if(hasLowerBound(qp.lower_x(i))) n++;
        }
    }
    return n;
}

void QPSwiftSolver::toQpSwift(const wbc::QuadraticProgram &qp){

    G.setZero();

    P = qp.H;
    c = qp.g;
    A = qp.A;
    b = qp.b;

    // create inequalities matrix (inequalities constraints + bounds), skipping unbounded sides
    uint row = 0;
    for(int i = 0; i < qp.nin; i++){
        if(hasUpperBound(qp.upper_y(i))){                                // Cx <= upper_y
            G.row(row) = qp.C.row(i);
            h(row++) = qp.upper_y(i);
        }
        if(hasLowerBound(qp.lower_y(i))){                                // -Cx <= -lower_y
            G.row(row) = -qp.C.row(i);
            h(row++) = -qp.lower_y(i);
        }
    }
    if(qp.bounded){
        for(int i = 0; i < qp.nq; i++){
            if(hasUpperBound(qp.upper_x(i))){                            // x <= upper_x
                G(row, i) = 1.0;
                h(row++) = qp.upper_x(i);
            }
            if(hasLowerBound(qp.lower_x(i))){                            // -x <= -lower_x
                G(row, i) = -1.0;
                h(row++) = -qp.lower_x(i);
            }
        }
    }
    assert(row == (uint)n_ineq);

    // QP_SETUP_dense allocates a new problem, so the one of the previous call has to be released
    if(my_qp)
        QP_CLEANUP_dense(my_qp);

    my_qp = QP_SETUP_dense(n_dec,                   // Number decision variables
                           n_ineq,                  // Number inequality constraints
                           n_eq,                    // Number equality constraints
                           (double*)P.data(),       // Hessian matrix
                           (double*)A.data(),       // Equality constraint matrix
                           (double*)G.data(),       // Inequality constraint matrix
                           (double*)c.data(),       // Cost function gradient vector
                           (double*)h.data(),       // Inequality constraint vector
                           (double*)b.data(),       // Equality constraint vector
                           NULL,
                           COLUMN_MAJOR_ORDERING);

    my_qp->options->maxit = options.maxit;
    my_qp->options->reltol = options.reltol;
    my_qp->options->abstol = options.abstol;
    my_qp->options->sigma = options.sigma;
    my_qp->options->verbose = options.verbose;
}

void QPSwiftSolver::solve(const wbc::HierarchicalQP &hierarchical_qp, Eigen::VectorXd &solver_output, bool allow_warm_start){

    assert(hierarchical_qp.size() == 1);
    const wbc::QuadraticProgram &qp = hierarchical_qp[0];
    assert(qp.isValid());

    if(!allow_warm_start)
        configured = false;

    // The number of inequality rows depends on how many sides are actually bounded, so it can
    // change between calls (e.g. when contacts are switched on or off) and the solver has to be
    // set up again when it does.
    const uint n_ineq_new = countInequalities(qp);
    if(n_dec != (int)qp.nq || n_eq != (int)qp.neq || n_ineq != (int)n_ineq_new)
        configured = false;

    if(!configured){
        // Count equality / inequality constraints
        n_dec = qp.nq;
        n_eq = qp.neq;
        n_ineq = n_ineq_new;

        A.resize(n_eq, n_dec);
        b.resize(n_eq);
        G.resize(n_ineq, n_dec);
        h.resize(n_ineq);
        P.resize(n_dec, n_dec);
        c.resize(n_dec);
        solver_output.resize(n_dec);
        configured = true;
    }

    toQpSwift(qp);

    qp_int exit_code = QP_SOLVE(my_qp);

    switch(exit_code){
    case QP_OPTIMAL:{
        break;
    }
    case QP_MAXIT:{
        throw std::runtime_error("QPSwiftSolver failed: Maximum Iterations reached");
    }
    case QP_FATAL:{
        throw std::runtime_error("QPSwiftSolver failed: Unknown error");
    }
    case QP_KKTFAIL:{
        throw std::runtime_error("QPSwiftSolver failed: LDL Factorization failed");
    }
    }

    for(int i = 0; i < n_dec; i++)
        solver_output[i] = my_qp->x[i];
}

}
