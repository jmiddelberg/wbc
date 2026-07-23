#include "ClarabelSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include <Eigen/Core>
#include <stdexcept>
#include <string>

namespace wbc {

QPSolverRegistry<ClarabelSolver> ClarabelSolver::reg("clarabel");

static std::string solverStatusToString(clarabel::SolverStatus status){
    switch(status){
        case clarabel::SolverStatus::Unsolved:               return "Unsolved";
        case clarabel::SolverStatus::Solved:                 return "Solved";
        case clarabel::SolverStatus::PrimalInfeasible:       return "PrimalInfeasible";
        case clarabel::SolverStatus::DualInfeasible:         return "DualInfeasible";
        case clarabel::SolverStatus::AlmostSolved:           return "AlmostSolved";
        case clarabel::SolverStatus::AlmostPrimalInfeasible: return "AlmostPrimalInfeasible";
        case clarabel::SolverStatus::AlmostDualInfeasible:   return "AlmostDualInfeasible";
        case clarabel::SolverStatus::MaxIterations:          return "MaxIterations";
        case clarabel::SolverStatus::MaxTime:                return "MaxTime";
        case clarabel::SolverStatus::NumericalError:         return "NumericalError";
        case clarabel::SolverStatus::InsufficientProgress:   return "InsufficientProgress";
        default:                                             return "Unknown";
    }
}

ClarabelSolver::ClarabelSolver() : _actual_n_iter(0)
{
    settings.verbose = false;
}

/// solve problem:
/// min  0.5 * x'Hx + g'x
/// s.t. Ax = b
///      lower_y <= Cx <= upper_y
///      lower_x <=  x <= upper_x   (only if qp.bounded)
///
/// transcribed to Clarabel's conic form  min 0.5 x'Px + q'x  s.t.  Âx + s = b̂, s in K:
///      [ A ] x            = b            -> ZeroCone(neq)
///      [ C ] x  <= upper_y                \
///      [-C ] x  <= -lower_y                > NonnegativeCone(2*nin + (bounded ? 2*nq : 0))
///      [ I ] x  <= upper_x (if bounded)   |
///      [-I ] x  <= -lower_x (if bounded)  /
void ClarabelSolver::solve(const wbc::HierarchicalQP& hierarchical_qp, Eigen::VectorXd& solver_output, bool allow_warm_start)
{
    (void)allow_warm_start; // Clarabel is an interior-point solver and is rebuilt on every call.

    assert(hierarchical_qp.size() == 1);

    const wbc::QuadraticProgram& qp = hierarchical_qp[0];
    assert(qp.isValid());

    const uint n_var = qp.nq;
    const uint n_bnd = qp.bounded ? qp.nq : 0;
    const uint m_zero = qp.neq;                         // rows in the zero (equality) cone
    const uint m_nn   = 2*qp.nin + 2*n_bnd;             // rows in the nonnegative (inequality) cone
    const uint m_total = m_zero + m_nn;

    // Objective: Clarabel reads the upper triangle of the (symmetric) Hessian.
    Eigen::MatrixXd H_upper = qp.H.triangularView<Eigen::Upper>();
    _P = H_upper.sparseView();
    _P.makeCompressed();

    _q = qp.g;

    // Stacked constraint matrix and right-hand side.
    _A_dense.resize(m_total, n_var);
    _A_dense.setZero();
    _b.resize(m_total);

    uint row = 0;
    if(qp.neq > 0){                                     // equalities  Ax = b
        _A_dense.middleRows(row, qp.neq) = qp.A;
        _b.segment(row, qp.neq) = qp.b;
        row += qp.neq;
    }
    if(qp.nin > 0){                                     // inequalities Cx <= upper_y  and  -Cx <= -lower_y
        _A_dense.middleRows(row, qp.nin) = qp.C;
        _b.segment(row, qp.nin) = qp.upper_y;
        row += qp.nin;
        _A_dense.middleRows(row, qp.nin) = -qp.C;
        _b.segment(row, qp.nin) = -qp.lower_y;
        row += qp.nin;
    }
    if(qp.bounded){                                     // bounds  x <= upper_x  and  -x <= -lower_x
        _A_dense.middleRows(row, qp.nq) = Eigen::MatrixXd::Identity(qp.nq, qp.nq);
        _b.segment(row, qp.nq) = qp.upper_x;
        row += qp.nq;
        _A_dense.middleRows(row, qp.nq) = -Eigen::MatrixXd::Identity(qp.nq, qp.nq);
        _b.segment(row, qp.nq) = -qp.lower_x;
        row += qp.nq;
    }

    _A = _A_dense.sparseView();
    _A.makeCompressed();

    std::vector<clarabel::SupportedConeT<double>> cones;
    if(m_zero > 0)
        cones.push_back(clarabel::ZeroConeT<double>(m_zero));
    if(m_nn > 0)
        cones.push_back(clarabel::NonnegativeConeT<double>(m_nn));

    clarabel::DefaultSolver<double> solver(_P, _q, _A, _b, cones, settings);
    solver.solve();

    clarabel::DefaultSolution<double> solution = solver.solution();

    solver_output.resize(n_var);
    solver_output = solution.x;

    _actual_n_iter = solution.iterations;

    if(solution.status != clarabel::SolverStatus::Solved &&
       solution.status != clarabel::SolverStatus::AlmostSolved){
        qp.print();
        throw std::runtime_error("Clarabel returned non-optimal status: " + solverStatusToString(solution.status));
    }
}

} // namespace wbc
