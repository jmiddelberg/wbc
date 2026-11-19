#include "AcadosOcpQPSolver.hpp"
#include "AcadosQPBackend.hpp"
#include "../../core/QuadraticProgram.hpp"

#include <acados/ocp_qp/ocp_qp_common.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp_sol.h>

#include <cstdlib>
#include <stdexcept>
#include <vector>

// The installed headers of these modules (acados/ocp_qp/ocp_qp_osqp.h, ocp_qp_clarabel.h) include the solver
// headers by their path in the acados source tree, which acados does not install. Only the entry point of each module
// is needed, so it is declared here. ACADOS_WITH_* come with the acados CMake target.
extern "C" {
#ifdef ACADOS_WITH_OSQP
void ocp_qp_osqp_config_initialize_default(void *config);
#endif
#ifdef ACADOS_WITH_CLARABEL
void ocp_qp_clarabel_config_initialize_default(void *config);
#endif
}

namespace wbc {

namespace {

/// The QP as a single stage acados OCP QP, solved directly by the OCP QP solver module
class OcpBackend : public acados_qp::Backend{
public:
    explicit OcpBackend(AcadosOcpQPSolver::Module module) : module(module){}
    ~OcpBackend(){ destroy(); }

    void create(const acados_qp::Dims& d) override{
        nv = d.nv;
        ne = d.ne;
        nb = d.nb;
        ng = d.ne + d.ng;   // general constraints: equalities, then inequalities

        config_mem = calloc(1, ocp_qp_solver_config_calculate_size());
        config = ocp_qp_solver_config_assign(config_mem);
        switch(module){
#ifdef ACADOS_WITH_OSQP
        case AcadosOcpQPSolver::MODULE_OSQP: ocp_qp_osqp_config_initialize_default(config); break;
#endif
#ifdef ACADOS_WITH_CLARABEL
        case AcadosOcpQPSolver::MODULE_CLARABEL: ocp_qp_clarabel_config_initialize_default(config); break;
#endif
        default: throw std::runtime_error("acados was built without the OCP QP solver module " + std::to_string(module));
        }

        dims_mem = calloc(1, ocp_qp_dims_calculate_size(0));
        dims = ocp_qp_dims_assign(0, dims_mem);
        d_ocp_qp_dim_set_nx(0, 0, dims);
        d_ocp_qp_dim_set_nu(0, nv, dims);
        d_ocp_qp_dim_set_nbx(0, 0, dims);
        d_ocp_qp_dim_set_nbu(0, nb, dims);
        d_ocp_qp_dim_set_ng(0, ng, dims);
        d_ocp_qp_dim_set_ns(0, 0, dims);
        d_ocp_qp_dim_set_nbxe(0, 0, dims);
        d_ocp_qp_dim_set_nbue(0, 0, dims);
        d_ocp_qp_dim_set_nge(0, ne, dims);

        qp_in_mem = calloc(1, ocp_qp_in_calculate_size(dims));
        qp_in = ocp_qp_in_assign(dims, qp_in_mem);
        qp_out_mem = calloc(1, ocp_qp_out_calculate_size(dims));
        qp_out = ocp_qp_out_assign(dims, qp_out_mem);

        opts_mem = calloc(1, config->opts_calculate_size(config, dims));
        opts = config->opts_assign(config, dims, opts_mem);
        config->opts_initialize_default(config, dims, opts);

        // The equalities are the first ne general constraints. The indices count the box constraints first.
        if(ne > 0){
            std::vector<int> idxe(ne);
            for(int i = 0; i < ne; i++)
                idxe[i] = nb + i;
            d_ocp_qp_set_idxe(0, idxe.data(), qp_in);
        }

        D.resize(ng, nv);
        lg.resize(ng);
        ug.resize(ng);
        lg_mask.resize(ng);
        ug_mask.resize(ng);
        lg_mask.head(ne).setOnes();
        ug_mask.head(ne).setOnes();
    }

    void createSolver() override{
        config->opts_update(config, dims, opts);
        mem_mem = calloc(1, config->memory_calculate_size(config, dims, opts));
        mem = config->memory_assign(config, dims, opts, mem_mem);
        const acados_size_t work_size = config->workspace_calculate_size(config, dims, opts);
        work = work_size > 0 ? calloc(1, work_size) : nullptr;
    }

    void destroy() override{
        // OSQP and Clarabel allocate their solver objects outside of the acados memory, terminate frees them
        if(mem && config && config->terminate)
            config->terminate(config, mem, work);
        for(void* p : {work, mem_mem, opts_mem, qp_out_mem, qp_in_mem, dims_mem, config_mem})
            free(p);
        work = mem_mem = mem = opts_mem = opts = qp_out_mem = qp_in_mem = dims_mem = config_mem = nullptr;
        qp_out = nullptr;
        qp_in = nullptr;
        dims = nullptr;
        config = nullptr;
    }

    bool created() const override{ return config != nullptr; }

    void setOption(const char* field, void* value) override{
        config->opts_set(config, opts, field, value);
    }

    void setData(const QuadraticProgram& qp, const acados_qp::Bounds& b) override{
        // All matrices are column-major, as Eigen stores them by default. HPIPM's setters only read the data, but
        // are not const-correct.
        d_ocp_qp_set_R(0, const_cast<double*>(qp.H.data()), qp_in);
        d_ocp_qp_set_r(0, const_cast<double*>(qp.g.data()), qp_in);

        if(nb > 0){
            d_ocp_qp_set_idxbu(0, const_cast<int*>(b.idxb.data()), qp_in);
            d_ocp_qp_set_lbu(0, const_cast<double*>(b.lb.data()), qp_in);
            d_ocp_qp_set_ubu(0, const_cast<double*>(b.ub.data()), qp_in);
            // Which sides are absent can change between two calls (contacts are switched on and off), so the masks
            // are written every time
            d_ocp_qp_set_lbu_mask(0, const_cast<double*>(b.lb_mask.data()), qp_in);
            d_ocp_qp_set_ubu_mask(0, const_cast<double*>(b.ub_mask.data()), qp_in);
        }

        if(ng > 0){
            const int nin = ng - ne;
            if(ne > 0){
                D.topRows(ne) = qp.A;
                lg.head(ne) = qp.b;
                ug.head(ne) = qp.b;
            }
            if(nin > 0){
                D.bottomRows(nin) = qp.C;
                lg.tail(nin) = Eigen::Map<const Eigen::VectorXd>(b.lg.data(), nin);
                ug.tail(nin) = Eigen::Map<const Eigen::VectorXd>(b.ug.data(), nin);
                lg_mask.tail(nin) = Eigen::Map<const Eigen::VectorXd>(b.lg_mask.data(), nin);
                ug_mask.tail(nin) = Eigen::Map<const Eigen::VectorXd>(b.ug_mask.data(), nin);
            }
            d_ocp_qp_set_D(0, D.data(), qp_in);
            d_ocp_qp_set_lg(0, lg.data(), qp_in);
            d_ocp_qp_set_ug(0, ug.data(), qp_in);
            d_ocp_qp_set_lg_mask(0, lg_mask.data(), qp_in);
            d_ocp_qp_set_ug_mask(0, ug_mask.data(), qp_in);
        }
    }

    int solve() override{
        return config->evaluate(config, qp_in, qp_out, opts, mem, work);
    }

    void getSolution(Eigen::VectorXd& x) const override{
        d_ocp_qp_sol_get_u(0, qp_out, x.data());
    }

    int iterations() const override{
        return qp_out ? static_cast<qp_info*>(qp_out->misc)->num_iter : 0;
    }

private:
    AcadosOcpQPSolver::Module module;
    int nv = 0, ne = 0, nb = 0, ng = 0;

    // acados keeps pointers between these objects, so they are allocated once per create() and do not move
    void *config_mem = nullptr, *dims_mem = nullptr, *qp_in_mem = nullptr, *qp_out_mem = nullptr;
    void *opts_mem = nullptr, *mem_mem = nullptr;
    qp_solver_config *config = nullptr;
    ocp_qp_dims *dims = nullptr;
    ocp_qp_in *qp_in = nullptr;
    ocp_qp_out *qp_out = nullptr;
    void *opts = nullptr;
    void *mem = nullptr;
    void *work = nullptr;

    // General constraints [A; C] and their bounds, the equalities with lg = ug = b
    Eigen::MatrixXd D;
    Eigen::VectorXd lg, ug, lg_mask, ug_mask;
};

}

AcadosOcpQPSolver::AcadosOcpQPSolver(Module module, const std::string& name, int default_warm_start, int max_warm_start) :
    AcadosQPSolver(std::unique_ptr<acados_qp::Backend>(new OcpBackend(module)), name, default_warm_start, max_warm_start){
}

}
