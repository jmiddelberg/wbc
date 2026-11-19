#ifndef WBC_SOLVERS_ACADOS_DENSE_QP_TEST_CASES_HPP
#define WBC_SOLVERS_ACADOS_DENSE_QP_TEST_CASES_HPP

// Test cases that every solver behind the acados dense QP interface has to pass. The test of each
// solver instantiates them with WBC_DENSE_QP_TEST_CASES(prefix, factory, tol), where factory
// creates a configured solver and tol is the accuracy that the solver reaches.

#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <functional>
#include <memory>
#include <stdexcept>
#include "core/QuadraticProgram.hpp"
#include "solvers/acados_dense_qp/AcadosDenseQPSolver.hpp"

namespace dense_qp_test{

using SolverFactory = std::function<std::unique_ptr<wbc::AcadosDenseQPSolver>()>;

// Task Jacobian used throughout the tests
inline Eigen::MatrixXd taskJacobian(){
    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    return A;
}

// Desired task space reference
inline Eigen::VectorXd taskReference(){
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;
    return y;
}

inline wbc::HierarchicalQP toHQP(const wbc::QuadraticProgram& qp){
    wbc::HierarchicalQP hqp;
    hqp.prios.push_back(qp);
    return hqp;
}

// Solve and report an exception as a test error, with the name of the solver
inline bool solve(wbc::AcadosDenseQPSolver& solver, const wbc::QuadraticProgram& qp, Eigen::VectorXd& x, bool allow_warm_start = true){
    BOOST_REQUIRE(qp.isValid());
    try{
        solver.solve(toHQP(qp), x, allow_warm_start);
    }
    catch(const std::exception& e){
        BOOST_ERROR(solver.getName() << ": " << e.what());
        return false;
    }
    BOOST_CHECK_EQUAL(solver.getReturnCode(), 0);
    BOOST_CHECK_EQUAL(x.size(), qp.nq);
    return true;
}

inline void checkFeasible(const wbc::QuadraticProgram& qp, const Eigen::VectorXd& x, double tol){
    if(qp.neq > 0)
        BOOST_CHECK_SMALL((qp.A*x - qp.b).cwiseAbs().maxCoeff(), tol);
    for(int j = 0; j < qp.nin; j++){
        const double c = (qp.C.row(j)*x).value();
        if(wbc::hasLowerBound(qp.lower_y(j)))
            BOOST_CHECK_MESSAGE(c >= qp.lower_y(j) - tol, "constraint " << j << ": " << c << " < " << qp.lower_y(j));
        if(wbc::hasUpperBound(qp.upper_y(j)))
            BOOST_CHECK_MESSAGE(c <= qp.upper_y(j) + tol, "constraint " << j << ": " << c << " > " << qp.upper_y(j));
    }
    if(qp.bounded){
        for(int j = 0; j < qp.nq; j++){
            if(wbc::hasLowerBound(qp.lower_x(j)))
                BOOST_CHECK_MESSAGE(x(j) >= qp.lower_x(j) - tol, "variable " << j << ": " << x(j) << " < " << qp.lower_x(j));
            if(wbc::hasUpperBound(qp.upper_x(j)))
                BOOST_CHECK_MESSAGE(x(j) <= qp.upper_x(j) + tol, "variable " << j << ": " << x(j) << " > " << qp.upper_x(j));
        }
    }
}

// min ||Ax-y||, no constraints: the task has to be solved exactly
inline void withoutConstraints(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 0, 0, false);
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x))
        BOOST_CHECK_SMALL((A*x - y).cwiseAbs().maxCoeff(), tol);
}

// min ||x|| s.t. Ax = y, with a square A: the solution is unique
inline void withEqualityConstraints(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 6, 0, false);
    qp.H.setIdentity();
    qp.g.setZero();
    qp.A = taskJacobian();
    qp.b = taskReference();

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x))
        checkFeasible(qp, x, tol);
}

// min ||x|| s.t. y-0.1 <= Ax <= y+0.1
inline void withInequalityConstraints(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 0, 6, false);
    qp.H.setIdentity();
    qp.g.setZero();
    qp.C = taskJacobian();
    qp.lower_y = taskReference() - Eigen::VectorXd::Constant(6, 1e-1);
    qp.upper_y = taskReference() + Eigen::VectorXd::Constant(6, 1e-1);

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x)){
        checkFeasible(qp, x, tol);
        // x = 0 is not feasible, so some constraints are active
        BOOST_CHECK(x.norm() > 1e-3);
    }
}

// min ||Ax-y|| with bounds that are far away: the task has to be solved exactly
inline void withInactiveBounds(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 0, 0, true);
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);
    qp.lower_x.setConstant(-1e10);
    qp.upper_x.setConstant(+1e10);

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x)){
        BOOST_CHECK_SMALL((A*x - y).cwiseAbs().maxCoeff(), tol);
        checkFeasible(qp, x, tol);
    }
}

// min ||Ax-y|| with bounds that become active
inline void withActiveBounds(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 0, 0, true);
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);
    qp.lower_x.setConstant(-0.1);
    qp.upper_x.setConstant(+0.1);

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x)){
        checkFeasible(qp, x, tol);
        BOOST_CHECK(x.cwiseAbs().maxCoeff() > 0.1 - 1e-6);
        BOOST_CHECK(solver->getIterations() > 0);
    }
}

// min ||x|| s.t. equalities, inactive inequalities and inactive bounds: the solution is the
// minimum norm solution of the equalities, which checks that each constraint type ends up in the
// right place
inline void withAllConstraintTypes(const SolverFactory& make, double tol){
    wbc::QuadraticProgram qp;
    qp.resize(6, 2, 2, true);
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.H.setIdentity();
    qp.g.setZero();
    qp.A = A.topRows(2);
    qp.b = y.head(2);
    qp.C = A.bottomRows(2);
    qp.lower_y.setConstant(-10);
    qp.upper_y.setConstant(+10);
    qp.lower_x.setConstant(-10);
    qp.upper_x.setConstant(+10);

    auto solver = make();
    Eigen::VectorXd x;
    if(solve(*solver, qp, x)){
        Eigen::VectorXd expected = qp.A.transpose()*(qp.A*qp.A.transpose()).inverse()*qp.b;
        BOOST_CHECK_SMALL((x - expected).cwiseAbs().maxCoeff(), tol);
        checkFeasible(qp, x, tol);
    }
}

// A sequence of QPs with the same dimensions, as in a control loop. Every warm started solution has
// to be equal to the solution of a new solver. Then a cold start of the same solver, and a change
// of the dimensions.
inline void repeatedSolve(const SolverFactory& make, double tol){
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    auto solver = make();
    Eigen::VectorXd x, x_ref;
    wbc::QuadraticProgram qp;

    for(int i = 0; i < 10; i++){
        // Moving reference and equality constraints, so that different constraints become active.
        // The equalities go through a point inside the bounds, so every QP is feasible.
        qp.resize(6, 2, 6, true);
        Eigen::VectorXd y_i = y + Eigen::VectorXd::Constant(6, 0.1*i);
        Eigen::VectorXd x_feasible = Eigen::VectorXd::Constant(6, 0.04*i);
        qp.H = A.transpose()*A + 1e-8*Eigen::MatrixXd::Identity(6,6);
        qp.g = -(A.transpose()*y_i);
        qp.A = A.topRows(2);
        qp.b = qp.A*x_feasible;
        qp.C = A;
        qp.lower_y.setConstant(-1.5);
        qp.upper_y.setConstant(+1.5);
        qp.lower_x.setConstant(-0.5);
        qp.upper_x.setConstant(+0.5);

        if(!solve(*solver, qp, x))
            return;
        checkFeasible(qp, x, tol);

        auto reference = make();
        if(solve(*reference, qp, x_ref))
            BOOST_CHECK_MESSAGE((x - x_ref).cwiseAbs().maxCoeff() < 10*tol, "cycle " << i << ": warm started solution differs by " << (x - x_ref).cwiseAbs().maxCoeff());
    }

    // A cold start of the same solver gives the same solution, and the next call warm starts again
    Eigen::VectorXd x_cold;
    if(solve(*solver, qp, x_cold, false))
        BOOST_CHECK_SMALL((x - x_cold).cwiseAbs().maxCoeff(), 10*tol);
    if(solve(*solver, qp, x_cold))
        BOOST_CHECK_SMALL((x - x_cold).cwiseAbs().maxCoeff(), 10*tol);

    // A change of the dimensions creates the solver again
    qp.resize(6, 0, 0, false);
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);
    if(solve(*solver, qp, x))
        BOOST_CHECK_SMALL((A*x - y).cwiseAbs().maxCoeff(), tol);
}

// Bounds marked with wbc::INF are absent. This must not change the solution: the same QP with
// bounds that are finite but too wide to become active has to give the same solution.
inline void unboundedConstraints(const SolverFactory& make, double tol){
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    // no_bound is either wbc::INF or a finite value far outside the feasible region
    auto makeQP = [&](double no_bound){
        wbc::QuadraticProgram qp;
        qp.resize(6, 2, 6, true);
        qp.H = A.transpose()*A;
        qp.H.diagonal().array() += 1e-6;
        qp.g = -(A.transpose()*y);
        qp.A = A.topRows(2);
        qp.b = y.head(2);
        qp.C = A;
        // one-sided inequalities: no lower bound, upper bound partly active
        qp.lower_y.setConstant(-no_bound);
        qp.upper_y = y + Eigen::VectorXd::Constant(6, 1e-1);
        // three free variables, three that are really bounded
        qp.lower_x.setConstant(-no_bound);
        qp.upper_x.setConstant(+no_bound);
        qp.lower_x.tail(3).setConstant(-0.4);
        qp.upper_x.tail(3).setConstant(+0.4);
        return qp;
    };

    wbc::QuadraticProgram qp_inf = makeQP(wbc::INF);
    wbc::QuadraticProgram qp_wide = makeQP(1e3);

    auto solver_inf = make();
    auto solver_wide = make();
    Eigen::VectorXd x_inf, x_wide;
    if(!solve(*solver_inf, qp_inf, x_inf) || !solve(*solver_wide, qp_wide, x_wide))
        return;

    BOOST_CHECK_SMALL((x_inf - x_wide).cwiseAbs().maxCoeff(), 1e-4);
    checkFeasible(qp_inf, x_inf, std::max(tol, 1e-6));
    // at least one of the real bounds is active, otherwise dropping every bound would pass too
    BOOST_CHECK(x_inf.tail(3).cwiseAbs().maxCoeff() > 0.4 - 1e-3);
}

// Which sides are absent changes between two calls with the same dimensions, as when contacts are
// switched on and off. The warm start must not keep a constraint that is absent now: every
// solution has to be equal to the solution of a new solver.
inline void changingMasks(const SolverFactory& make, double tol){
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    // The constraints are built around a feasible point, so every QP is feasible: masking a side
    // only removes constraints. The cost pulls the solution away from it, so some constraints and
    // bounds are active.
    const Eigen::VectorXd x_feasible = Eigen::VectorXd::Constant(6, 0.1);

    auto solver = make();
    Eigen::VectorXd x, x_ref;

    for(int i = 0; i < 8; i++){
        wbc::QuadraticProgram qp;
        qp.resize(6, 1, 6, true);
        qp.H = A.transpose()*A + 1e-6*Eigen::MatrixXd::Identity(6,6);
        qp.g = -(A.transpose()*(2*y));
        qp.A = A.topRows(1);
        qp.b = qp.A*x_feasible;
        qp.C = A;
        // even cycles: all bounds and upper constraint limits, odd cycles: half of them absent
        const bool even = (i % 2 == 0);
        const Eigen::VectorXd upper_y = A*x_feasible + Eigen::VectorXd::Constant(6, 0.05);
        for(int j = 0; j < 6; j++){
            const bool absent = !even && (j % 2 == (i/2) % 2);
            qp.lower_x(j) = absent ? -wbc::INF : -0.3;
            qp.upper_x(j) = absent ?  wbc::INF :  0.3;
            qp.lower_y(j) = -wbc::INF;
            qp.upper_y(j) = absent ?  wbc::INF :  upper_y(j);
        }

        if(!solve(*solver, qp, x))
            return;
        checkFeasible(qp, x, std::max(tol, 1e-6));

        auto reference = make();
        if(solve(*reference, qp, x_ref))
            BOOST_CHECK_MESSAGE((x - x_ref).cwiseAbs().maxCoeff() < 1e-5, "cycle " << i << ": warm started solution differs by " << (x - x_ref).cwiseAbs().maxCoeff());
    }
}

// An infeasible QP throws and gives a return code other than ACADOS_SUCCESS. The next feasible QP
// with the same dimensions is solved again, so a failed solve must not spoil the warm start.
inline void failureAndRecovery(const SolverFactory& make, double tol){
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    wbc::QuadraticProgram qp;
    qp.resize(6, 1, 0, true);
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);
    qp.A.setZero();
    qp.A(0,0) = 1.0;
    qp.lower_x.setConstant(-0.5);
    qp.upper_x.setConstant(+0.5);

    auto solver = make();
    solver->setMaxIter(1000);
    Eigen::VectorXd x;

    qp.b(0) = 0.2;  // feasible
    solve(*solver, qp, x);

    qp.b(0) = 2.0;  // infeasible: x0 = 2 and x0 <= 0.5
    BOOST_CHECK_THROW(solver->solve(toHQP(qp), x), std::runtime_error);
    BOOST_CHECK(solver->getReturnCode() != 0);

    qp.b(0) = 0.2;  // feasible again
    if(solve(*solver, qp, x))
        checkFeasible(qp, x, std::max(tol, 1e-6));
}

// Option validation that is common to all solvers
inline void commonOptions(const SolverFactory& make, int default_warm_start, int max_warm_start){
    auto solver = make();
    BOOST_CHECK_EQUAL(solver->getWarmStart(), default_warm_start);
    BOOST_CHECK_THROW(solver->setWarmStart(-1), std::invalid_argument);
    BOOST_CHECK_THROW(solver->setWarmStart(max_warm_start+1), std::invalid_argument);
    for(int level = 0; level <= max_warm_start; level++)
        BOOST_CHECK_NO_THROW(solver->setWarmStart(level));
    BOOST_CHECK_THROW(solver->setMaxIter(0), std::invalid_argument);
    BOOST_CHECK(solver->setTolerance(1e-6));
    BOOST_CHECK_NO_THROW(solver->setMaxIter(50));
}

// A tighter tolerance gives a solution that is at least as accurate. The active bounds problem is
// solved with a tight tolerance as reference.
inline void tolerance(const SolverFactory& make){
    wbc::QuadraticProgram qp;
    qp.resize(6, 0, 0, true);
    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y);
    qp.lower_x.setConstant(-0.4);
    qp.upper_x.setConstant(+0.4);

    auto reference = make();
    reference->setTolerance(1e-10);
    Eigen::VectorXd x_ref;
    if(!solve(*reference, qp, x_ref))
        return;

    for(double tol : {1e-3, 1e-6, 1e-9}){
        auto solver = make();
        BOOST_CHECK(solver->setTolerance(tol));
        Eigen::VectorXd x;
        if(!solve(*solver, qp, x))
            continue;
        checkFeasible(qp, x, tol);
        // A residual tolerance bounds the error of the solution only up to the conditioning of the
        // problem, hence the factor
        BOOST_CHECK_MESSAGE((x - x_ref).cwiseAbs().maxCoeff() < std::max(1e3*tol, 1e-6),
                            "tol " << tol << ": error " << (x - x_ref).cwiseAbs().maxCoeff());
    }
}

} // namespace dense_qp_test

/// Instantiate all common test cases for one solver
#define WBC_DENSE_QP_TEST_CASES(prefix, factory, tol) \
    BOOST_AUTO_TEST_CASE(prefix##_without_constraints){ dense_qp_test::withoutConstraints(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_with_equality_constraints){ dense_qp_test::withEqualityConstraints(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_with_inequality_constraints){ dense_qp_test::withInequalityConstraints(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_with_inactive_bounds){ dense_qp_test::withInactiveBounds(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_with_active_bounds){ dense_qp_test::withActiveBounds(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_with_all_constraint_types){ dense_qp_test::withAllConstraintTypes(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_repeated_solve){ dense_qp_test::repeatedSolve(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_unbounded_constraints){ dense_qp_test::unboundedConstraints(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_changing_masks){ dense_qp_test::changingMasks(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_failure_and_recovery){ dense_qp_test::failureAndRecovery(factory, tol); } \
    BOOST_AUTO_TEST_CASE(prefix##_tolerance){ dense_qp_test::tolerance(factory); }

#endif
