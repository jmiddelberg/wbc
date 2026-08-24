#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sys/time.h>
#include <Eigen/Dense>
#include "../../../core/QuadraticProgram.hpp"
#include "../DAQPSolver.hpp"

using namespace wbc;
using namespace std;

// DAQP terminates in a KKT point of the (regularized) problem, so the solution is exact up to
// the conditioning of the problem
static const double TOL = 1e-9;

// Task Jacobian used throughout the tests
static Eigen::MatrixXd taskJacobian(){
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
static Eigen::VectorXd taskReference(){
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;
    return y;
}

BOOST_AUTO_TEST_CASE(solver_daqp_without_constraints)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = false;

    // Solve the problem min(||Ax-b||) without constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    Eigen::VectorXd test = A*solver_output;
    for(uint j = 0; j < NO_JOINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);
}

BOOST_AUTO_TEST_CASE(solver_daqp_with_equality_constraints)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 6;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = false;

    // Solve the problem min(||x||), subject Ax=b --> encode the task as constraint
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = I  and g = 0

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    qp.g.setZero();
    qp.H.setIdentity();

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.A = A;
    qp.b = y;

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    Eigen::VectorXd test = A*solver_output;
    for(uint j = 0; j < NO_EQ_CONSTRAINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);
}

BOOST_AUTO_TEST_CASE(solver_daqp_with_inequality_constraints)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 6;
    const bool WITH_BOUNDS = false;

    // Solve the problem min(||x||), subject to lb <= Ax <= ub --> encode the task as inequality constraint
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = I  and g = 0

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    qp.g.setZero();
    qp.H.setIdentity();

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();
    qp.C = A;
    qp.lower_y = y - Eigen::VectorXd::Constant(qp.nq, 1e-1);
    qp.upper_y = y + Eigen::VectorXd::Constant(qp.nq, 1e-1);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    Eigen::VectorXd test = A*solver_output;
    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp.lower_y(j)-TOL) <= test(j) && test(j) <= (qp.upper_y(j)+TOL));
}

BOOST_AUTO_TEST_CASE(solver_daqp_bounded)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = true;

    // Solve the problem min(||Ax-b||) with bound constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    qp.lower_x.setConstant(-1e10);
    qp.upper_x.setConstant(+1e10);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    // The bounds are inactive here, so the task has to be solved exactly
    Eigen::VectorXd test = A*solver_output;
    for(uint j = 0; j < NO_JOINTS; ++j){
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);
        BOOST_CHECK((qp.lower_x(j)-TOL) <= solver_output(j) && solver_output(j) <= (qp.upper_x(j)+TOL));
    }
}

BOOST_AUTO_TEST_CASE(solver_daqp_with_active_bounds)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = true;

    // Same problem as above, but with bounds that are tight enough to become active

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    qp.lower_x.setConstant(-0.1);
    qp.upper_x.setConstant(+0.1);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp.lower_x(j)-TOL) <= solver_output(j) && solver_output(j) <= (qp.upper_x(j)+TOL));

    // At least one bound has to be active, otherwise the test case would be pointless
    BOOST_CHECK(solver_output.cwiseAbs().maxCoeff() > 0.1 - TOL);
}

BOOST_AUTO_TEST_CASE(solver_daqp_with_all_constraint_types)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 2;
    const int NO_IN_CONSTRAINTS = 2;
    const bool WITH_BOUNDS = true;

    // Solve the problem min(||x||) subject to equality constraints, inequality constraints and bounds
    // at the same time. This checks that all constraint types end up in the right place in the
    // stacked problem that is passed to DAQP.

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    qp.g.setZero();
    qp.H.setIdentity();

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    qp.A = A.topRows(NO_EQ_CONSTRAINTS);
    qp.b = y.head(NO_EQ_CONSTRAINTS);
    qp.C = A.bottomRows(NO_IN_CONSTRAINTS);
    // Inequalities and bounds are chosen wide enough to stay inactive
    qp.lower_y.setConstant(-10);
    qp.upper_y.setConstant(+10);
    qp.lower_x.setConstant(-10);
    qp.upper_x.setConstant(+10);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    DAQPSolver solver;

    Eigen::VectorXd solver_output;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    // With only the equalities active, the solution is the minimum norm solution of A*x = b
    Eigen::MatrixXd A_eq = qp.A;
    Eigen::VectorXd expected = A_eq.transpose()*(A_eq*A_eq.transpose()).inverse()*qp.b;
    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK(fabs(solver_output(j) - expected(j)) < TOL);

    Eigen::VectorXd test = qp.C*solver_output;
    for(uint j = 0; j < NO_IN_CONSTRAINTS; ++j)
        BOOST_CHECK((qp.lower_y(j)-TOL) <= test(j) && test(j) <= (qp.upper_y(j)+TOL));
}

BOOST_AUTO_TEST_CASE(solver_daqp_repeated_solve)
{
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 2;
    const int NO_IN_CONSTRAINTS = 6;
    const bool WITH_BOUNDS = true;

    // Solve a sequence of problems with the same dimensions, which is what happens in a control
    // loop: the solver re-uses its workspace and warm starts from the previous active set.

    Eigen::MatrixXd A = taskJacobian();
    Eigen::VectorXd y = taskReference();

    DAQPSolver solver;
    Eigen::VectorXd solver_output;

    wbc::QuadraticProgram qp;
    wbc::HierarchicalQP hqp;

    for(uint i = 0; i < 10; i++){
        qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

        // Moving task reference and moving equality constraint, so that different constraints
        // become active over time. The equalities are constructed from a point inside the bounds,
        // which keeps the problem feasible in every cycle.
        Eigen::VectorXd y_i = y + Eigen::VectorXd::Constant(NO_JOINTS, 0.1*i);
        Eigen::VectorXd x_feasible = Eigen::VectorXd::Constant(NO_JOINTS, 0.04*i);

        qp.H = A.transpose()*A + 1e-8*Eigen::MatrixXd::Identity(NO_JOINTS,NO_JOINTS);
        qp.g = -(A.transpose()*y_i).transpose();
        qp.A = A.topRows(NO_EQ_CONSTRAINTS);
        qp.b = qp.A*x_feasible;
        qp.C = A.bottomRows(NO_IN_CONSTRAINTS);
        qp.lower_y.setConstant(-1.5);
        qp.upper_y.setConstant(+1.5);
        qp.lower_x.setConstant(-0.5);
        qp.upper_x.setConstant(+0.5);

        BOOST_CHECK(qp.isValid());
        hqp.prios.clear();
        hqp << qp;

        BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
        BOOST_CHECK(solver.getExitFlag() == DAQP_EXIT_OPTIMAL);

        // Check that the solution of every single cycle is feasible
        Eigen::VectorXd test_eq = qp.A*solver_output;
        for(uint j = 0; j < NO_EQ_CONSTRAINTS; ++j)
            BOOST_CHECK(fabs(test_eq(j) - qp.b(j)) < TOL);

        Eigen::VectorXd test_in = qp.C*solver_output;
        for(uint j = 0; j < NO_IN_CONSTRAINTS; ++j)
            BOOST_CHECK((qp.lower_y(j)-TOL) <= test_in(j) && test_in(j) <= (qp.upper_y(j)+TOL));

        for(uint j = 0; j < NO_JOINTS; ++j)
            BOOST_CHECK((qp.lower_x(j)-TOL) <= solver_output(j) && solver_output(j) <= (qp.upper_x(j)+TOL));
    }

    // Cold starting the last problem has to give the same solution as the warm started solve
    Eigen::VectorXd solver_output_cold;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output_cold, false));
    BOOST_CHECK((solver_output - solver_output_cold).norm() < TOL);

    // A change of the problem dimensions has to trigger a re-configuration of the solver
    qp.resize(NO_JOINTS, 0, 0, false);
    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    BOOST_CHECK(qp.isValid());
    hqp.prios.clear();
    hqp << qp;

    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));

    Eigen::VectorXd test = A*solver_output;
    for(uint j = 0; j < NO_JOINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);
}

BOOST_AUTO_TEST_CASE(solver_daqp_unbounded_constraints)
{
    // Bounds marked with wbc::INF are absent. The wrapper has to hand that to the solver in the way
    // the solver expects, and doing so must not change the solution: solving the same problem with
    // bounds that are finite but far too wide to ever become active has to give the same answer.
    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 2;
    const int NO_IN_CONSTRAINTS = 6;

    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;

    // no_bound is either wbc::INF or a finite value that is far outside the feasible region
    auto makeQP = [&](double no_bound){
        wbc::QuadraticProgram qp;
        qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, true);
        qp.H = A.transpose()*A;
        qp.H.diagonal().array() += 1e-6;
        qp.g = -(A.transpose()*y).transpose();
        qp.A = A.topRows(NO_EQ_CONSTRAINTS);
        qp.b = y.head(NO_EQ_CONSTRAINTS);
        qp.C = A;
        // one-sided inequalities: no lower bound, upper bound partly active
        qp.lower_y.setConstant(-no_bound);
        qp.upper_y = y + Eigen::VectorXd::Constant(NO_IN_CONSTRAINTS, 1e-1);
        // three free variables, three that are really bounded
        qp.lower_x.setConstant(-no_bound);
        qp.upper_x.setConstant(+no_bound);
        qp.lower_x.tail(3).setConstant(-0.4);
        qp.upper_x.tail(3).setConstant(+0.4);
        return qp;
    };

    wbc::QuadraticProgram qp_inf = makeQP(wbc::INF);
    wbc::QuadraticProgram qp_wide = makeQP(1e3);
    BOOST_CHECK(qp_inf.isValid());

    wbc::HierarchicalQP hqp_inf, hqp_wide;
    hqp_inf << qp_inf;
    hqp_wide << qp_wide;

    Eigen::VectorXd out_inf, out_wide;
    DAQPSolver solver_inf, solver_wide;
    BOOST_CHECK_NO_THROW(solver_inf.solve(hqp_inf, out_inf));
    BOOST_CHECK_NO_THROW(solver_wide.solve(hqp_wide, out_wide));

    BOOST_CHECK_EQUAL(out_inf.size(), NO_JOINTS);
    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK_SMALL(out_inf(j) - out_wide(j), 1e-4);

    // the bounds that are not marked as absent must still be enforced ...
    for(uint j = 3; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp_inf.lower_x(j)-1e-4) <= out_inf(j) && out_inf(j) <= (qp_inf.upper_x(j)+1e-4));
    // ... and at least one of them has to be active, otherwise dropping every bound would pass too
    BOOST_CHECK(out_inf.tail(3).cwiseAbs().maxCoeff() > 0.4 - 1e-3);

    // the equality constraints are never affected by the substitution
    BOOST_CHECK_SMALL((qp_inf.A*out_inf - qp_inf.b).cwiseAbs().maxCoeff(), 1e-4);
}
