#include "AcadosDenseQPSolver.hpp"
#include "AcadosQPBackend.hpp"
#include "../../core/QuadraticProgram.hpp"

#include <acados_c/dense_qp_interface.h>
#include <acados/dense_qp/dense_qp_common.h>
#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_sol.h>

#include <cstdlib>

namespace wbc {

namespace {

/// The QP as an acados dense QP, solved with dense_qp_solve()
class DenseBackend : public acados_qp::Backend{
public:
    explicit DenseBackend(int acados_solver){
        plan.qp_solver = static_cast<dense_qp_solver_t>(acados_solver);
        dims.nv = dims.ne = dims.nb = dims.ng = dims.ns = -1;
    }
    ~DenseBackend(){ destroy(); }

    void create(const acados_qp::Dims& d) override{
        dims.nv = d.nv;
        dims.ne = d.ne;
        dims.nb = d.nb;
        dims.ng = d.ng;
        dims.ns = 0;
        config = dense_qp_config_create(&plan);
        qp_in = dense_qp_in_create(config, &dims);
        qp_out = dense_qp_out_create(config, &dims);
        opts = dense_qp_opts_create(config, &dims);
    }

    void createSolver() override{
        solver = dense_qp_create(config, &dims, opts);
    }

    void destroy() override{
        if(solver){
            // ProxQP allocates its solver object outside of the acados memory, terminate frees it
            if(config && config->terminate)
                config->terminate(config, solver->mem, solver->work);
            free(solver);
        }
        free(qp_out);
        free(qp_in);
        free(opts);
        free(config);
        solver = nullptr;
        qp_out = nullptr;
        qp_in = nullptr;
        opts = nullptr;
        config = nullptr;
    }

    bool created() const override{ return config != nullptr; }

    void setOption(const char* field, void* value) override{
        config->opts_set(config, opts, field, value);
    }

    void setData(const QuadraticProgram& qp, const acados_qp::Bounds& b) override{
        // All matrices are column-major, as Eigen stores them by default. HPIPM's setters only read the data, but
        // are not const-correct.
        d_dense_qp_set_all(const_cast<double*>(qp.H.data()),
                           const_cast<double*>(qp.g.data()),
                           const_cast<double*>(qp.A.data()),
                           const_cast<double*>(qp.b.data()),
                           const_cast<int*>(b.idxb.data()),
                           const_cast<double*>(b.lb.data()),
                           const_cast<double*>(b.ub.data()),
                           const_cast<double*>(qp.C.data()),
                           const_cast<double*>(b.lg.data()),
                           const_cast<double*>(b.ug.data()),
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,  // no slacks
                           qp_in);

        // d_dense_qp_set_all does not touch the masks. Which sides are absent can change between two calls
        // (contacts are switched on and off), so they are written every time.
        if(dims.nb > 0){
            d_dense_qp_set_lb_mask(const_cast<double*>(b.lb_mask.data()), qp_in);
            d_dense_qp_set_ub_mask(const_cast<double*>(b.ub_mask.data()), qp_in);
        }
        if(dims.ng > 0){
            d_dense_qp_set_lg_mask(const_cast<double*>(b.lg_mask.data()), qp_in);
            d_dense_qp_set_ug_mask(const_cast<double*>(b.ug_mask.data()), qp_in);
        }
    }

    int solve() override{
        return dense_qp_solve(solver, qp_in, qp_out);
    }

    void getSolution(Eigen::VectorXd& x) const override{
        d_dense_qp_sol_get_v(qp_out, x.data());
    }

    int iterations() const override{
        return qp_out ? static_cast<qp_info*>(qp_out->misc)->num_iter : 0;
    }

private:
    dense_qp_solver_plan plan;
    dense_qp_dims dims;             // acados keeps a pointer to this struct, so the backend must not move
    qp_solver_config *config = nullptr;
    dense_qp_in *qp_in = nullptr;
    dense_qp_out *qp_out = nullptr;
    void *opts = nullptr;
    dense_qp_solver *solver = nullptr;
};

}

AcadosDenseQPSolver::AcadosDenseQPSolver(int acados_solver, const std::string& name, int default_warm_start, int max_warm_start) :
    AcadosQPSolver(std::unique_ptr<acados_qp::Backend>(new DenseBackend(acados_solver)), name, default_warm_start, max_warm_start){
}

}
