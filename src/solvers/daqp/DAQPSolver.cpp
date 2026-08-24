#include "DAQPSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include <daqp/utils.h>
#include <stdexcept>
#include <string>

namespace wbc {

QPSolverRegistry<DAQPSolver> DAQPSolver::reg("daqp");

// The problem data changes completely between two calls to solve(), only the dimensions stay the same
static const int update_mask = DAQP_UPDATE_Rinv | DAQP_UPDATE_v | DAQP_UPDATE_M | DAQP_UPDATE_d | DAQP_UPDATE_sense;

static std::string exitFlagToString(int exitflag){
    switch(exitflag){
        case DAQP_EXIT_SOFT_OPTIMAL:            return "SoftOptimal";
        case DAQP_EXIT_OPTIMAL:                 return "Optimal";
        case DAQP_EXIT_INFEASIBLE:              return "Infeasible";
        case DAQP_EXIT_CYCLE:                   return "Cycle";
        case DAQP_EXIT_UNBOUNDED:               return "Unbounded";
        case DAQP_EXIT_ITERLIMIT:               return "IterationLimit";
        case DAQP_EXIT_NONCONVEX:               return "NonConvex";
        case DAQP_EXIT_OVERDETERMINED_INITIAL:  return "OverdeterminedInitialWorkingSet";
        case DAQP_EXIT_TIMELIMIT:               return "TimeLimit";
#ifdef DAQP_EXIT_UNSUPPORTED   // not defined by DAQP versions before v0.9
        case DAQP_EXIT_UNSUPPORTED:             return "Unsupported";
#endif
        default:                                return "Unknown (" + std::to_string(exitflag) + ")";
    }
}

DAQPSolver::DAQPSolver() : _work(), _result(), _problem(), _n_var(-1), _n_cstr(-1), _n_bnd(-1)
{
    daqp_default_settings(&_settings);
}

DAQPSolver::~DAQPSolver(){
    freeWorkspace();
}

void DAQPSolver::freeWorkspace(){
    // DAQP frees the settings that the workspace points to, so hide the ones owned by this class
    _work.settings = 0;
    free_daqp_workspace(&_work);
    free_daqp_ldp(&_work);
    _n_var = _n_cstr = _n_bnd = -1;
}

/// solve problem:
/// min  0.5 * x'Hx + g'x
/// s.t. Ax = b
///      lower_y <= Cx <= upper_y
///      lower_x <=  x <= upper_x   (only if qp.bounded)
///
/// transcribed to the form  min 0.5 x'Hx + f'x  s.t.  blower <= [I; A_daqp] x <= bupper  expected by
/// DAQP, where the first ms rows are the simple bounds on the first ms variables:
///      lower_x <=   x  <= upper_x  (if bounded, ms = nq rows)
///            b <=  Ax  <= b        (equalities, marked active and immutable)
///      lower_y <=  Cx  <= upper_y
void DAQPSolver::solve(const wbc::HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start)
{
    assert(hierarchical_qp.size() == 1);

    const wbc::QuadraticProgram& qp = hierarchical_qp[0];
    assert(qp.isValid());

    const int n_var  = qp.nq;
    const int n_bnd  = qp.bounded ? qp.nq : 0;   // simple bounds, DAQP expects them on the leading variables
    const int n_gen  = qp.neq + qp.nin;          // rows of the general constraint matrix
    const int n_cstr = n_bnd + n_gen;

    if(!allow_warm_start)
        configured = false;
    // The workspace is allocated for a fixed problem size, so it has to be rebuilt if that changes
    if(n_var != _n_var || n_cstr != _n_cstr || n_bnd != _n_bnd)
        configured = false;

    // Constraint matrix and bounds. DAQP reads the constraint matrix row-wise and expects the bounds
    // of the simple bounds before those of the general constraints.
    _A.resize(n_gen, n_var);
    _A.topRows(qp.neq) = qp.A;
    _A.bottomRows(qp.nin) = qp.C;

    // DAQP has its own value for an unbounded constraint, so wbc::INF is mapped onto it. A bound
    // beyond DAQP_INF is never a candidate for the working set, whereas a large finite one is.
    auto toDaqpInf = [](const Eigen::VectorXd& v){
        return v.unaryExpr([](double b){
            if(b >= INF)  return  (double)DAQP_INF;
            if(b <= -INF) return -(double)DAQP_INF;
            return b;
        }).eval();
    };

    _b_upper.resize(n_cstr);
    _b_lower.resize(n_cstr);
    _b_upper.head(n_bnd) = toDaqpInf(qp.upper_x);
    _b_lower.head(n_bnd) = toDaqpInf(qp.lower_x);
    _b_upper.segment(n_bnd, qp.neq) = qp.b;
    _b_lower.segment(n_bnd, qp.neq) = qp.b;
    _b_upper.tail(qp.nin) = toDaqpInf(qp.upper_y);
    _b_lower.tail(qp.nin) = toDaqpInf(qp.lower_y);

    // All constraints are inequalities, apart from the equalities, which DAQP encodes as constraints
    // that are active from the start and may never leave the working set.
    _sense.setZero(n_cstr);
    if(configured){
        // Warm start: initialize the working set with the constraints that were active in the
        // previous solution (the sign of the dual solution tells at which bound they were active)
        for(int i = 0; i < n_cstr; i++){
            if(_lam(i) > 0.0)
                _sense(i) = DAQP_ACTIVE;
            else if(_lam(i) < 0.0)
                _sense(i) = DAQP_ACTIVE | DAQP_LOWER;
        }
    }
    _sense.segment(n_bnd, qp.neq).setConstant(DAQP_ACTIVE | DAQP_IMMUTABLE);

    // DAQP does not modify the problem data, but its interface is not const-correct. The Hessian is
    // read row-wise, which makes no difference here since it is symmetric.
    _problem.n = n_var;
    _problem.m = n_cstr;
    _problem.ms = n_bnd;
    _problem.H = const_cast<double*>(qp.H.data());
    _problem.f = const_cast<double*>(qp.g.data());
    _problem.A = n_gen > 0 ? _A.data() : 0;
    _problem.bupper = _b_upper.data();
    _problem.blower = _b_lower.data();
    _problem.sense = _sense.data();

    if(!configured){
        _x.resize(n_var);
        _lam.setZero(n_cstr);
        _result.x = _x.data();
        _result.lam = _lam.data();

        freeWorkspace();
        _work.settings = &_settings;   // DAQP allocates its own settings if this is not set
        int flag = setup_daqp(&_problem, &_work, 0);
        if(flag < 0){
            qp.print();
            throw std::runtime_error("DAQP setup failed with exit flag: " + exitFlagToString(flag));
        }
        _n_var = n_var;
        _n_cstr = n_cstr;
        _n_bnd = n_bnd;
        configured = true;
    }
    else{
        int flag = daqp_update_ldp(update_mask, &_work, &_problem);
        if(flag < 0){
            qp.print();
            throw std::runtime_error("DAQP update failed with exit flag: " + exitFlagToString(flag));
        }
    }

    daqp_solve(&_result, &_work);

    if(_result.exitflag < 0){
        qp.print();
        throw std::runtime_error("DAQP returned non-optimal exit flag: " + exitFlagToString(_result.exitflag));
    }

    solver_output = _x;
}

} // namespace wbc
