#include "EiquadprogSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include <Eigen/Core>
#include <iostream>

namespace wbc {

QPSolverRegistry<EiquadprogSolver> EiquadprogSolver::reg("eiquadprog");

EiquadprogSolver::EiquadprogSolver()
{
    _n_iter = 100;
}

EiquadprogSolver::~EiquadprogSolver()
{

}

void EiquadprogSolver::solve(const wbc::HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool /*allow_warm_start*/)
{

    assert(hierarchical_qp.size() == 1);

    const wbc::QuadraticProgram &qp = hierarchical_qp[0];
    assert(qp.isValid());

    // Eiquadprog expects one-sided inequalities CI*x + ci0 >= 0, so every two-sided constraint has
    // to be split into two rows. Sides marked with wbc::INF are not emitted: eiquadprog has no
    // notion of an infinite bound and would propagate it into its dual variables, and even a large
    // finite bound is not free, since every row enters the CI*x product of every iteration.
    size_t n_in = 0;
    for(int i = 0; i < qp.nin; i++){
        if(hasLowerBound(qp.lower_y(i))) n_in++;
        if(hasUpperBound(qp.upper_y(i))) n_in++;
    }
    if(qp.bounded){
        for(int i = 0; i < qp.nq; i++){
            if(hasLowerBound(qp.lower_x(i))) n_in++;
            if(hasUpperBound(qp.upper_x(i))) n_in++;
        }
    }
    size_t n_eq = qp.neq;
    size_t n_var = qp.nq;

    // The number of rows depends on how many sides are actually bounded, so it can change between
    // calls (e.g. when contacts are switched on or off) and the workspace has to be rebuilt
    if(_n_in != n_in || _n_eq != n_eq || _n_var != n_var)
        configured = false;

    if(!configured) 
    {
        _solver.reset(n_var, n_eq, n_in);
        _solver.setMaxIter(_n_iter);

        // hessian and gradient are ok (don#t need to be stacked)
        // equality contraint is ok also
        // configuring  inequalities constraints matrices
        _CI_mtx.resize(n_in, n_var);
        _ci0_vec.resize(n_in);

        _n_in = n_in;
        _n_eq = n_eq;
        _n_var = n_var;
        configured = true;
    }

    _CI_mtx.setZero();

    // create inequalities constraint matrix (inequalities + bounds), skipping unbounded sides
    size_t row = 0;
    for(int i = 0; i < qp.nin; i++){
        if(hasLowerBound(qp.lower_y(i))){                    // Cx - lower_y >= 0
            _CI_mtx.row(row) = qp.C.row(i);
            _ci0_vec(row++) = -qp.lower_y(i);
        }
        if(hasUpperBound(qp.upper_y(i))){                    // -Cx + upper_y >= 0
            _CI_mtx.row(row) = -qp.C.row(i);
            _ci0_vec(row++) = qp.upper_y(i);
        }
    }
    if(qp.bounded){
        for(int i = 0; i < qp.nq; i++){
            if(hasLowerBound(qp.lower_x(i))){                // x - lower_x >= 0
                _CI_mtx(row, i) = 1.0;
                _ci0_vec(row++) = -qp.lower_x(i);
            }
            if(hasUpperBound(qp.upper_x(i))){                // -x + upper_x >= 0
                _CI_mtx(row, i) = -1.0;
                _ci0_vec(row++) = qp.upper_x(i);
            }
        }
    }
    assert(row == n_in);

    namespace eq = eiquadprog::solvers;

    Eigen::VectorXd out(qp.nq);

    eq::EiquadprogFast_status status = _solver.solve_quadprog(
        qp.H, qp.g, qp.A, -qp.b, _CI_mtx, _ci0_vec, out);
    
    solver_output.resize(qp.nq);
    solver_output = out;

    if(status == eq::EiquadprogFast_status::EIQUADPROG_FAST_UNBOUNDED){
        //qp.print();
        throw std::runtime_error("Eiquadprog returned error status:unbounded.");
    }
    if(status == eq::EiquadprogFast_status::EIQUADPROG_FAST_MAX_ITER_REACHED){
        //qp.print();
        throw std::runtime_error("Eiquadprog returned error status: max iterations reached.");
    }
    if(status == eq::EiquadprogFast_status::EIQUADPROG_FAST_REDUNDANT_EQUALITIES){
        //qp.print();
        throw std::runtime_error("Eiquadprog returned error status: redundant equalities.");
    }
    if(status == eq::EiquadprogFast_status::EIQUADPROG_FAST_INFEASIBLE){
        //qp.print();
        throw std::runtime_error("Eiquadprog returned error status: infeasible.");
    }

    _actual_n_iter = _solver.getIteratios();
}
}
