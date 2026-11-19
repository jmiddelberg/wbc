#include "AcadosDenseQPSolver.hpp"
#include "../../core/QuadraticProgram.hpp"

#include <acados_c/dense_qp_interface.h>
#include <acados/dense_qp/dense_qp_common.h>
#include <acados/utils/types.h>
#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_sol.h>

#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace wbc {

/// Everything that has an acados type, so that the header does not need the acados include paths
struct AcadosDenseQPSolver::Acados{
    dense_qp_solver_plan plan;
    dense_qp_dims dims;             // acados keeps a pointer to this struct, so it must not move
    qp_solver_config *config = nullptr;
    dense_qp_in *qp_in = nullptr;
    dense_qp_out *qp_out = nullptr;
    void *opts = nullptr;
    dense_qp_solver *solver = nullptr;
    int applied_warm_start = -1;    // value of the acados option warm_start in opts

    // Bounds in the layout acados expects. A side that is absent is masked (mask 0.0); its value is
    // never read by the solvers, but is kept finite.
    std::vector<int> idxb;
    std::vector<double> lb, ub, lb_mask, ub_mask;
    std::vector<double> lg, ug, lg_mask, ug_mask;
};

AcadosDenseQPSolver::AcadosDenseQPSolver(int acados_solver, const std::string& name, int default_warm_start, int max_warm_start) :
    acados(new Acados()),
    name(name),
    max_warm_start(max_warm_start),
    warm_start(default_warm_start),
    tolerance(-1.0),
    iter_max(-1),
    return_code(ACADOS_READY){
    acados->plan.qp_solver = static_cast<dense_qp_solver_t>(acados_solver);
    acados->dims.nv = acados->dims.ne = acados->dims.nb = acados->dims.ng = acados->dims.ns = -1;
    configured = false;
}

AcadosDenseQPSolver::~AcadosDenseQPSolver(){
    destroy();
}

void AcadosDenseQPSolver::destroy(){
    Acados &a = *acados;
    if(a.solver){
        // ProxQP allocates its solver object outside of the acados memory, terminate frees it
        if(a.config && a.config->terminate)
            a.config->terminate(a.config, a.solver->mem, a.solver->work);
        free(a.solver);
    }
    free(a.qp_out);
    free(a.qp_in);
    free(a.opts);
    free(a.config);
    a.solver = nullptr;
    a.qp_out = nullptr;
    a.qp_in = nullptr;
    a.opts = nullptr;
    a.config = nullptr;
    a.applied_warm_start = -1;
}

void AcadosDenseQPSolver::setOption(const std::string& field, void* value){
    Acados &a = *acados;
    if(!a.config || !a.opts)
        throw std::logic_error(name + ": setOption() is only valid inside setSolverOptions()");
    a.config->opts_set(a.config, a.opts, field.c_str(), value);
}

void AcadosDenseQPSolver::create(){
    Acados &a = *acados;

    a.config = dense_qp_config_create(&a.plan);
    a.qp_in = dense_qp_in_create(a.config, &a.dims);
    a.qp_out = dense_qp_out_create(a.config, &a.dims);
    a.opts = dense_qp_opts_create(a.config, &a.dims);

    // The same fields and iteration limit as the MPC (model_predictive_control.cpp), so that both
    // are comparable
    if(tolerance > 0){
        for(const char *field : {"tol_stat", "tol_eq", "tol_ineq", "tol_comp", "tol_dual_gap"})
            setOption(field, &tolerance);
        int iter_max_tol = 1000;
        setOption("iter_max", &iter_max_tol);
    }
    if(iter_max > 0)
        setOption("iter_max", &iter_max);
    setOption("warm_start", &warm_start);
    a.applied_warm_start = warm_start;

    setSolverOptions();

    // Some solvers size their memory from the options, so the memory is created after them
    a.solver = dense_qp_create(a.config, &a.dims, a.opts);
}

/// The masks are what tells the solvers that a side is absent: every acados dense QP solver
/// replaces a masked bound by its own infinity before it uses it. The value of a masked side is
/// set to +-wbc::INF, which keeps it finite.
void AcadosDenseQPSolver::setData(const QuadraticProgram& qp){
    Acados &a = *acados;
    const int nb = a.dims.nb;
    const int ng = a.dims.ng;

    for(int i = 0; i < nb; i++){
        const bool has_lb = hasLowerBound(qp.lower_x(i));
        const bool has_ub = hasUpperBound(qp.upper_x(i));
        a.lb[i] = has_lb ? qp.lower_x(i) : -INF;
        a.ub[i] = has_ub ? qp.upper_x(i) : INF;
        a.lb_mask[i] = has_lb ? 1.0 : 0.0;
        a.ub_mask[i] = has_ub ? 1.0 : 0.0;
    }
    for(int i = 0; i < ng; i++){
        const bool has_lg = hasLowerBound(qp.lower_y(i));
        const bool has_ug = hasUpperBound(qp.upper_y(i));
        a.lg[i] = has_lg ? qp.lower_y(i) : -INF;
        a.ug[i] = has_ug ? qp.upper_y(i) : INF;
        a.lg_mask[i] = has_lg ? 1.0 : 0.0;
        a.ug_mask[i] = has_ug ? 1.0 : 0.0;
    }

    // All matrices are column-major, as Eigen stores them by default. HPIPM's setters only read
    // the data, but are not const-correct.
    d_dense_qp_set_all(const_cast<double*>(qp.H.data()),
                       const_cast<double*>(qp.g.data()),
                       const_cast<double*>(qp.A.data()),
                       const_cast<double*>(qp.b.data()),
                       a.idxb.data(),
                       a.lb.data(),
                       a.ub.data(),
                       const_cast<double*>(qp.C.data()),
                       a.lg.data(),
                       a.ug.data(),
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  // no slacks
                       a.qp_in);

    // d_dense_qp_set_all does not touch the masks. Which sides are absent can change between two
    // calls (contacts are switched on and off), so they are written every time.
    if(nb > 0){
        d_dense_qp_set_lb_mask(a.lb_mask.data(), a.qp_in);
        d_dense_qp_set_ub_mask(a.ub_mask.data(), a.qp_in);
    }
    if(ng > 0){
        d_dense_qp_set_lg_mask(a.lg_mask.data(), a.qp_in);
        d_dense_qp_set_ug_mask(a.ug_mask.data(), a.qp_in);
    }
}

void AcadosDenseQPSolver::solve(const HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start){

    assert(hierarchical_qp.size() == 1);
    const QuadraticProgram &qp = hierarchical_qp[0];
    assert(qp.isValid());

    Acados &a = *acados;

    const int nv = qp.nq;
    const int ne = qp.neq;
    const int nb = qp.bounded ? qp.nq : 0;
    const int ng = qp.nin;

    // The acados memory is sized for fixed dimensions
    if(nv != a.dims.nv || ne != a.dims.ne || nb != a.dims.nb || ng != a.dims.ng)
        configured = false;

    if(!configured){
        destroy();
        a.dims.nv = nv;
        a.dims.ne = ne;
        a.dims.nb = nb;
        a.dims.ng = ng;
        a.dims.ns = 0;
        a.idxb.resize(nb);
        for(int i = 0; i < nb; i++)
            a.idxb[i] = i;
        a.lb.resize(nb);
        a.ub.resize(nb);
        a.lb_mask.resize(nb);
        a.ub_mask.resize(nb);
        a.lg.resize(ng);
        a.ug.resize(ng);
        a.lg_mask.resize(ng);
        a.ug_mask.resize(ng);
        create();
        configured = true;
    }

    // A cold start only for this call, the memory and the options stay as they are
    int ws = allow_warm_start ? warm_start : 0;
    if(ws != a.applied_warm_start){
        a.config->opts_set(a.config, a.opts, "warm_start", &ws);
        a.applied_warm_start = ws;
    }

    setData(qp);

    return_code = dense_qp_solve(a.solver, a.qp_in, a.qp_out);
    if(return_code != ACADOS_SUCCESS)
        throw std::runtime_error(name + " returned: " + returnCodeToString(return_code));

    solver_output.resize(nv);
    d_dense_qp_sol_get_v(a.qp_out, solver_output.data());
}

bool AcadosDenseQPSolver::setTolerance(double tol){
    tolerance = tol;
    reset();
    return true;
}

void AcadosDenseQPSolver::setMaxIter(int n){
    if(n < 1)
        throw std::invalid_argument(name + ": the iteration limit must be at least 1, got " + std::to_string(n));
    iter_max = n;
    reset();
}

/// Validated here, because acados terminates the process on a warm start level it does not know,
/// and some valid levels assume constant matrices, which a WBC QP never has
void AcadosDenseQPSolver::setWarmStart(int level){
    if(level < 0 || level > max_warm_start)
        throw std::invalid_argument(name + ": invalid warm start level " + std::to_string(level) +
                                    ", must be between 0 and " + std::to_string(max_warm_start));
    warm_start = level;
}

int AcadosDenseQPSolver::getIterations() const{
    if(!acados->qp_out)
        return 0;
    return static_cast<qp_info*>(acados->qp_out->misc)->num_iter;
}

/// The acados return values (return_values_t), not the status of the solver: acados translates it
std::string AcadosDenseQPSolver::returnCodeToString(int code){
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
