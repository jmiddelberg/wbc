#include "HPIPMSolver.hpp"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include <hpipm_d_dense_qp.h>
#include <iostream>

namespace wbc{

QPSolverRegistry<HPIPMSolver> HPIPMSolver::reg("hpipm");

HPIPMSolver::HPIPMSolver(){
    plan.qp_solver = DENSE_QP_HPIPM;
    qp_in = 0;
    opts = 0;
    qp_out = 0;
    qp_solver = 0;
    config = 0;
}

HPIPMSolver::~HPIPMSolver(){
    if(qp_in)
        free(qp_in);
    if(opts)
        free(opts);
    if(qp_out)
        free(qp_out);
    if(qp_solver)
        free(qp_solver);
    if(config)
        free(config );
}

/// HPIPM does not need a bound for every variable: dims.nb says how many variables are bounded and
/// idxb says which ones. It also carries a mask over the constraint sides, where a 0.0 switches that
/// side off. Both matter here because HPIPM is an interior-point method - a bound that is present
/// but far away is not free, it contributes a slack, a barrier term and a complementarity condition
/// to every iteration.
///
/// Variables that are unbounded on both sides are therefore left out of idxb entirely, variables
/// that are bounded on one side only are included with the absent side masked off, and the same is
/// done for the one-sided rows of the general constraints (the friction cone facets, for example).
bool HPIPMSolver::updateBounds(const wbc::QuadraticProgram &qp){

    // Only the size matters for the return value: d_dense_qp_set_all rewrites the contents of idxb
    // on every call, so a different set of the same size needs no reconfiguration
    const size_t nb_before = idxb.size();

    idxb.clear();
    lb.clear();
    ub.clear();
    lb_mask.clear();
    ub_mask.clear();

    if(qp.bounded){
        for(int i = 0; i < qp.nq; i++){
            const bool has_lb = hasLowerBound(qp.lower_x(i));
            const bool has_ub = hasUpperBound(qp.upper_x(i));
            if(!has_lb && !has_ub)      // free variable, HPIPM never has to know about it
                continue;
            idxb.push_back(i);
            // The value of a masked side is not used, but it must stay finite: HPIPM still reads it
            lb.push_back(has_lb ? qp.lower_x(i) : 0.0);
            ub.push_back(has_ub ? qp.upper_x(i) : 0.0);
            lb_mask.push_back(has_lb ? 1.0 : 0.0);
            ub_mask.push_back(has_ub ? 1.0 : 0.0);
        }
    }

    lg_mask.resize(qp.nin);
    ug_mask.resize(qp.nin);
    for(int i = 0; i < qp.nin; i++){
        lg_mask[i] = hasLowerBound(qp.lower_y(i)) ? 1.0 : 0.0;
        ug_mask[i] = hasUpperBound(qp.upper_y(i)) ? 1.0 : 0.0;
    }

    return idxb.size() != nb_before;
}

void HPIPMSolver::solve(const HierarchicalQP &hierarchical_qp, Eigen::VectorXd &solver_output, bool allow_warm_start){

    assert(hierarchical_qp.size() == 1);
    const QuadraticProgram &qp = hierarchical_qp[0];
    assert(qp.isValid());

    if(!allow_warm_start)
        configured = false;

    // The number of bounded variables can change between calls, and it is part of the dimensions
    // HPIPM allocates for, so the problem has to be created again when it does
    if(updateBounds(qp))
        configured = false;

    if(!configured){

        dims.nv = qp.nq;
        dims.ne = qp.neq;
        dims.nb = idxb.size();
        dims.ng = qp.nin;
        dims.ns = 0;

        if(config)
            free(config);
        if(qp_in)
            free(qp_in);
        if(opts)
            free(opts);

        config = dense_qp_config_create(&plan);
        qp_in = dense_qp_in_create(config, &dims);
        opts = dense_qp_opts_create(config, &dims);

        dense_qp_hpipm_opts *hpipm_opts = (dense_qp_hpipm_opts *)opts;
        hpipm_opts->hpipm_opts->warm_start = 0;
        hpipm_opts->hpipm_opts->mode = SPEED;

        if(qp_out)
            free(qp_out);
        if(qp_solver)
            free(qp_solver);

        qp_out = dense_qp_out_create(config, &dims);
        qp_solver = dense_qp_create(config, &dims, opts);

        configured = true;
    }
    else{
        dense_qp_hpipm_opts *hpipm_opts = (dense_qp_hpipm_opts *)opts;
        hpipm_opts->hpipm_opts->warm_start = 1;
    }

    d_dense_qp_set_all((double* )qp.H.data(),
                       (double* )qp.g.data(),
                       (double* )qp.A.data(),
                       (double* )qp.b.data(),
                       idxb.data(),
                       lb.data(),
                       ub.data(),
                       (double* )qp.C.data(),
                       (double* )qp.lower_y.data(),
                       (double* )qp.upper_y.data(),
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       qp_in);

    // d_dense_qp_set_all does not touch the mask, so it is applied afterwards. The mask persists
    // across calls, but which sides are absent can change (contacts are switched on and off), so it
    // is rewritten every time.
    if(dims.nb > 0){
        d_dense_qp_set_lb_mask(lb_mask.data(), qp_in);
        d_dense_qp_set_ub_mask(ub_mask.data(), qp_in);
    }
    if(dims.ng > 0){
        d_dense_qp_set_lg_mask(lg_mask.data(), qp_in);
        d_dense_qp_set_ug_mask(ug_mask.data(), qp_in);
    }

    int ret = dense_qp_solve(qp_solver, qp_in, qp_out);
    if(ret != ACADOS_SUCCESS)
        throw std::runtime_error("HPIPM returned: " + returnCodeToString(ret));

    solver_output.resize(qp.nq);
    solver_output.setZero();
    d_dense_qp_sol_get_v(qp_out,solver_output.data());
}

void HPIPMSolver::setOptions(std::string &field,  void *value){
    dense_qp_hpipm_opts *hpipm_opts = (dense_qp_hpipm_opts *)opts;
    d_dense_qp_ipm_arg_set((char*)field.c_str(), value, hpipm_opts->hpipm_opts);
}

std::string HPIPMSolver::returnCodeToString(int code){
    switch(code){
    case SUCCESS: return "Found solution satisfying accuracy tolerance";
    case MAX_ITER: return "Maximum iteration number reached";
    case MIN_STEP: return "Minimum step length reached";
    case NAN_SOL: return "NaN in solution detected";
    case INCONS_EQ: return "unconsistent equality constraints";
    default: return "Unknown error code";
    }
}

}
