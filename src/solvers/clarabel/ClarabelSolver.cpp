#include "ClarabelSolver.hpp"
#include "../../core/QuadraticProgram.hpp"
#include <Eigen/Core>
#include <limits>
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

// Clarabel.rs pins the values of the linear solver enum that crosses the FFI boundary
// (AUTO = 0, QDLDL = 1, FAER = 2, MKL = 3, PANUA = 4), while the enum in
// clarabel/cpp/DefaultSettings.hpp only declares the methods that were compiled in and numbers
// them consecutively. Both agree as long as the cargo features are enabled in that same order,
// which is what the assertions below check: a build in which they disagree (e.g. pardiso-panua
// without faer-sparse) fails to compile instead of silently selecting the wrong method.
#ifdef FEATURE_FAER_SPARSE
static_assert(static_cast<int>(clarabel::ClarabelDirectSolveMethods::FAER) == 2,
              "clarabel::ClarabelDirectSolveMethods::FAER does not have the value expected by Clarabel.rs");
#endif
#ifdef FEATURE_PARDISO_MKL
static_assert(static_cast<int>(clarabel::ClarabelDirectSolveMethods::PARDISO_MKL) == 3,
              "clarabel::ClarabelDirectSolveMethods::PARDISO_MKL does not have the value expected by Clarabel.rs. "
              "Clarabel has to be built with the faer-sparse feature as well");
#endif
#ifdef FEATURE_PARDISO_PANUA
static_assert(static_cast<int>(clarabel::ClarabelDirectSolveMethods::PARDISO_PANUA) == 4,
              "clarabel::ClarabelDirectSolveMethods::PARDISO_PANUA does not have the value expected by Clarabel.rs. "
              "Clarabel has to be built with the faer-sparse and pardiso-mkl features as well");
#endif

static clarabel::ClarabelDirectSolveMethods toClarabelMethod(ClarabelLinearSolver method){
    switch(method){
        case clarabel_linsolver_auto:            return clarabel::ClarabelDirectSolveMethods::AUTO;
        case clarabel_linsolver_qdldl:           return clarabel::ClarabelDirectSolveMethods::QDLDL;
#ifdef FEATURE_FAER_SPARSE
        case clarabel_linsolver_faer:            return clarabel::ClarabelDirectSolveMethods::FAER;
#endif
#ifdef FEATURE_PARDISO_MKL
        case clarabel_linsolver_pardiso_mkl:     return clarabel::ClarabelDirectSolveMethods::PARDISO_MKL;
#endif
#ifdef FEATURE_PARDISO_PANUA
        case clarabel_linsolver_pardiso_panua:   return clarabel::ClarabelDirectSolveMethods::PARDISO_PANUA;
#endif
        default: break;
    }
    throw std::runtime_error("Clarabel linear system solver '" + ClarabelSolver::linearSolverName(method) +
                             "' is not available in this build");
}

static ClarabelLinearSolver fromClarabelMethod(clarabel::ClarabelDirectSolveMethods method){
    switch(method){
        case clarabel::ClarabelDirectSolveMethods::AUTO:            return clarabel_linsolver_auto;
        case clarabel::ClarabelDirectSolveMethods::QDLDL:           return clarabel_linsolver_qdldl;
#ifdef FEATURE_FAER_SPARSE
        case clarabel::ClarabelDirectSolveMethods::FAER:            return clarabel_linsolver_faer;
#endif
#ifdef FEATURE_PARDISO_MKL
        case clarabel::ClarabelDirectSolveMethods::PARDISO_MKL:     return clarabel_linsolver_pardiso_mkl;
#endif
#ifdef FEATURE_PARDISO_PANUA
        case clarabel::ClarabelDirectSolveMethods::PARDISO_PANUA:   return clarabel_linsolver_pardiso_panua;
#endif
    }
    throw std::runtime_error("Clarabel reported the unknown linear system solver " +
                             std::to_string(static_cast<int>(method)));
}

// Cargo feature and wbc cmake option that a linear system solver requires, used to tell the user
// what to do when the requested method is not part of this build.
static std::string linearSolverBuildHint(ClarabelLinearSolver method){
    switch(method){
        case clarabel_linsolver_faer:
            return "the 'faer-sparse' cargo feature and -DCLARABEL_FEATURE_FAER_SPARSE=ON";
        case clarabel_linsolver_pardiso_mkl:
            return "the 'pardiso-mkl' cargo feature and -DCLARABEL_FEATURE_PARDISO_MKL=ON";
        case clarabel_linsolver_pardiso_panua:
            return "the 'pardiso-panua' cargo feature and -DCLARABEL_FEATURE_PARDISO_PANUA=ON";
        default:
            return "";
    }
}

ClarabelSolver::ClarabelSolver() : _actual_n_iter(0)
{
    settings.verbose = false;
    settings.direct_solve_method = toClarabelMethod(linear_solver);
}

void ClarabelSolver::setOptions(clarabel::DefaultSettings<double> opt){
    settings = opt;
    // Keep getLinearSolver() in sync with the options that are actually in effect
    linear_solver = fromClarabelMethod(opt.direct_solve_method);
}

void ClarabelSolver::setLinearSolver(ClarabelLinearSolver method){
    if(!linearSolverAvailable(method)){
        std::string available;
        for(const std::string& name : availableLinearSolvers())
            available += (available.empty() ? "" : ", ") + name;
        throw std::runtime_error("Clarabel linear system solver '" + linearSolverName(method) + "' is not available "
                                 "in this build (available: " + available + "). Rebuild libclarabel_c with " +
                                 linearSolverBuildHint(method) + ", see "
                                 "https://clarabel.org/stable/user_guide_linsolvers/");
    }
    linear_solver = method;
    settings.direct_solve_method = toClarabelMethod(method);
}

void ClarabelSolver::setLinearSolver(const std::string& method){
    setLinearSolver(linearSolverFromName(method));
}

bool ClarabelSolver::linearSolverAvailable(ClarabelLinearSolver method){
    switch(method){
        case clarabel_linsolver_auto:
        case clarabel_linsolver_qdldl:
            return true;
        case clarabel_linsolver_faer:
#ifdef FEATURE_FAER_SPARSE
            return true;
#else
            return false;
#endif
        case clarabel_linsolver_pardiso_mkl:
#ifdef FEATURE_PARDISO_MKL
            return true;
#else
            return false;
#endif
        case clarabel_linsolver_pardiso_panua:
#ifdef FEATURE_PARDISO_PANUA
            return true;
#else
            return false;
#endif
    }
    return false;
}

std::vector<std::string> ClarabelSolver::availableLinearSolvers(){
    const ClarabelLinearSolver all[] = {clarabel_linsolver_auto,
                                        clarabel_linsolver_qdldl,
                                        clarabel_linsolver_faer,
                                        clarabel_linsolver_pardiso_mkl,
                                        clarabel_linsolver_pardiso_panua};
    std::vector<std::string> names;
    for(ClarabelLinearSolver method : all){
        if(linearSolverAvailable(method))
            names.push_back(linearSolverName(method));
    }
    return names;
}

std::string ClarabelSolver::linearSolverName(ClarabelLinearSolver method){
    switch(method){
        case clarabel_linsolver_auto:          return "auto";
        case clarabel_linsolver_qdldl:         return "qdldl";
        case clarabel_linsolver_faer:          return "faer";
        case clarabel_linsolver_pardiso_mkl:   return "pardiso-mkl";
        case clarabel_linsolver_pardiso_panua: return "pardiso-panua";
    }
    throw std::invalid_argument("Invalid Clarabel linear system solver " +
                                std::to_string(static_cast<int>(method)));
}

ClarabelLinearSolver ClarabelSolver::linearSolverFromName(const std::string& name){
    if(name == "auto")          return clarabel_linsolver_auto;
    if(name == "qdldl")         return clarabel_linsolver_qdldl;
    if(name == "faer")          return clarabel_linsolver_faer;
    if(name == "pardiso-mkl")   return clarabel_linsolver_pardiso_mkl;
    if(name == "pardiso-panua") return clarabel_linsolver_pardiso_panua;
    throw std::invalid_argument("Unknown Clarabel linear system solver '" + name + "'. Valid names are "
                                "auto, qdldl, faer, pardiso-mkl and pardiso-panua");
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
///
/// Right-hand sides of the NonnegativeCone block that reach wbc::INF, i.e. bounds the scenes marked
/// as absent, are replaced by an actual infinity so that Clarabel's presolve removes those rows
/// before solving (see setInfinity()).
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

    // Promote the bounds the scenes marked as absent (wbc::INF) to actual infinities. Clarabel's
    // presolve then drops these rows, which typically removes more than half of the
    // NonnegativeCone block and roughly halves the iteration count. Only the NonnegativeCone block
    // is touched, the equality rows are left alone. The substitution requires
    // settings.presolve_enable (which is on by default): with presolve disabled the infinities
    // would be fed to the solver directly.
    if(m_nn > 0 && settings.presolve_enable){
        _b.tail(m_nn) = (_b.tail(m_nn).array() >= infinity)
                            .select(std::numeric_limits<double>::infinity(), _b.tail(m_nn));
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

    // Which linear system solver clarabel_linsolver_auto resolved to
    _linsolver_used = fromClarabelMethod(solver.info().linsolver.name);

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
