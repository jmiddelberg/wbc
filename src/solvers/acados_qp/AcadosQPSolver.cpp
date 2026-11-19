#include "AcadosQPSolver.hpp"
#include "AcadosQPBackend.hpp"
#include "../../core/QuadraticProgram.hpp"

#include <acados/utils/types.h>

#include <stdexcept>

namespace wbc {

namespace acados_qp {

void Bounds::resize(const Dims& dims){
    idxb.resize(dims.nb);
    for(int i = 0; i < dims.nb; i++)
        idxb[i] = i;
    for(std::vector<double>* v : {&lb, &ub, &lb_mask, &ub_mask})
        v->resize(dims.nb);
    for(std::vector<double>* v : {&lg, &ug, &lg_mask, &ug_mask})
        v->resize(dims.ng);
}

/// The masks are what tells the solvers that a side is absent: every acados solver replaces a masked bound by its
/// own way to drop it before it uses it. The value of a masked side is set to +-wbc::INF, which keeps it finite.
void Bounds::set(const QuadraticProgram& qp){
    for(size_t i = 0; i < idxb.size(); i++){
        const bool has_lb = hasLowerBound(qp.lower_x(i));
        const bool has_ub = hasUpperBound(qp.upper_x(i));
        lb[i] = has_lb ? qp.lower_x(i) : -INF;
        ub[i] = has_ub ? qp.upper_x(i) : INF;
        lb_mask[i] = has_lb ? 1.0 : 0.0;
        ub_mask[i] = has_ub ? 1.0 : 0.0;
    }
    for(size_t i = 0; i < lg.size(); i++){
        const bool has_lg = hasLowerBound(qp.lower_y(i));
        const bool has_ug = hasUpperBound(qp.upper_y(i));
        lg[i] = has_lg ? qp.lower_y(i) : -INF;
        ug[i] = has_ug ? qp.upper_y(i) : INF;
        lg_mask[i] = has_lg ? 1.0 : 0.0;
        ug_mask[i] = has_ug ? 1.0 : 0.0;
    }
}

}

/// Everything that depends on the dimensions of the QP
struct AcadosQPSolver::State{
    acados_qp::Dims dims;
    acados_qp::Bounds bounds;
};

AcadosQPSolver::AcadosQPSolver(std::unique_ptr<acados_qp::Backend> backend, const std::string& name, int default_warm_start, int max_warm_start) :
    backend(std::move(backend)),
    state(new State()),
    name(name),
    max_warm_start(max_warm_start),
    warm_start(default_warm_start),
    applied_warm_start(-1),
    tolerance(-1.0),
    iter_max(-1),
    return_code(ACADOS_READY),
    setting_options(false){
    configured = false;
}

AcadosQPSolver::~AcadosQPSolver(){
    backend->destroy();
}

void AcadosQPSolver::setOption(const std::string& field, void* value){
    if(!setting_options)
        throw std::logic_error(name + ": setOption() is only valid inside setSolverOptions()");
    backend->setOption(field.c_str(), value);
}

void AcadosQPSolver::create(){
    backend->destroy();
    applied_warm_start = -1;
    state->bounds.resize(state->dims);
    backend->create(state->dims);

    setting_options = true;
    try{
        // The same fields and iteration limit as the MPC (model_predictive_control.cpp), so that both are comparable
        if(tolerance > 0){
            for(const char *field : {"tol_stat", "tol_eq", "tol_ineq", "tol_comp", "tol_dual_gap"})
                setOption(field, &tolerance);
            int iter_max_tol = 1000;
            setOption("iter_max", &iter_max_tol);
        }
        if(iter_max > 0)
            setOption("iter_max", &iter_max);
        setOption("warm_start", &warm_start);
        applied_warm_start = warm_start;

        setSolverOptions();
    }
    catch(...){
        setting_options = false;
        throw;
    }
    setting_options = false;

    backend->createSolver();
}

void AcadosQPSolver::solve(const HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start){

    assert(hierarchical_qp.size() == 1);
    const QuadraticProgram &qp = hierarchical_qp[0];
    assert(qp.isValid());

    acados_qp::Dims dims;
    dims.nv = qp.nq;
    dims.ne = qp.neq;
    dims.nb = qp.bounded ? qp.nq : 0;
    dims.ng = qp.nin;

    // The acados memory is sized for fixed dimensions
    if(dims != state->dims)
        configured = false;

    if(!configured){
        state->dims = dims;
        try{
            create();
        }
        catch(...){
            // invalid dimensions, so that the next call creates the memory again
            state->dims = acados_qp::Dims();
            throw;
        }
        configured = true;
    }

    // A cold start only for this call, the memory and the options stay as they are
    int ws = allow_warm_start ? warm_start : 0;
    if(ws != applied_warm_start){
        backend->setOption("warm_start", &ws);
        applied_warm_start = ws;
    }

    state->bounds.set(qp);
    backend->setData(qp, state->bounds);

    return_code = backend->solve();
    if(return_code != ACADOS_SUCCESS)
        throw std::runtime_error(name + " returned: " + returnCodeToString(return_code));

    solver_output.resize(dims.nv);
    backend->getSolution(solver_output);
}

bool AcadosQPSolver::setTolerance(double tol){
    tolerance = tol;
    reset();
    return true;
}

void AcadosQPSolver::setMaxIter(int n){
    if(n < 1)
        throw std::invalid_argument(name + ": the iteration limit must be at least 1, got " + std::to_string(n));
    iter_max = n;
    reset();
}

/// Validated here, because acados terminates the process on a warm start level it does not know,
/// and some valid levels assume constant matrices, which a WBC QP never has
void AcadosQPSolver::setWarmStart(int level){
    if(level < 0 || level > max_warm_start)
        throw std::invalid_argument(name + ": invalid warm start level " + std::to_string(level) +
                                    ", must be between 0 and " + std::to_string(max_warm_start));
    warm_start = level;
}

int AcadosQPSolver::getIterations() const{
    if(!backend->created())
        return 0;
    return backend->iterations();
}

/// The acados return values (return_values_t), not the status of the solver: acados translates it
std::string AcadosQPSolver::returnCodeToString(int code){
    switch(code){
    case ACADOS_SUCCESS: return "Found solution satisfying accuracy tolerance";
    case ACADOS_NAN_DETECTED: return "NaN in QP data or solution detected";
    case ACADOS_MAXITER: return "Maximum iteration number reached";
    case ACADOS_MINSTEP: return "Minimum step length reached";
    case ACADOS_QP_FAILURE: return "QP solver failed";
    case ACADOS_UNBOUNDED: return "Problem is unbounded";
    case ACADOS_TIMEOUT: return "Time limit reached";
    case ACADOS_INFEASIBLE: return "Problem is infeasible";
    case ACADOS_UNKNOWN: return "Unknown solver status";
    default: return "Solver specific status " + std::to_string(code);
    }
}

}
