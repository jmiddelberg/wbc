#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <sys/time.h>
#include "../../../core/QuadraticProgram.hpp"
#include "../EiquadprogSolver.hpp"

using namespace wbc;
using namespace std;

BOOST_AUTO_TEST_CASE(solver_qp_oases_without_constraints)
{
    srand (time(NULL));

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = false;
    const int NO_WSR = 20;

    // Solve the problem min(||Ax-b||) without constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T
    // For a 6x6 Constraint matrix this is approx. 3-5 times faster than encoding the task as constraint as below
    // With warm start, this solver is much faster (approx. 5 times) than in the initial run

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    qp.lower_x.resize(0);
    qp.upper_x.resize(0);
    qp.lower_y.resize(0);
    qp.upper_y.resize(0);
    qp.A.setIdentity();
    // Task Jacobian
    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    // Desired task space reference
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;

    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    EiquadprogSolver solver;
    solver.setMaxNIter(NO_WSR);

    BOOST_CHECK(solver.getMaxNIter() == NO_WSR);

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);
    //long useconds = end.tv_usec - start.tv_usec;

    /*cout<<"\n----------------------- Test Results ----------------------"<<endl<<endl;
    std::cout<<"Solver took "<<useconds<<" us "<<std::endl;
    cout<<"No of joints: "<<NO_JOINTS<<endl;
    cout<<"No of constraints: "<<NO_CONSTRAINTS<<endl;

    cout<<"\nSolver Input:"<<endl;
    cout<<"Constraint Matrix A:"<<endl; cout<<A<<endl;
    cout<<"Reference: y = "<<y.transpose()<<endl;

    cout<<"\nSolver Output: q_dot = "<<solver_output.transpose()<<endl;*/
    Eigen::VectorXd test = A*solver_output;
    //cout<<"Test: A * q_dot = "<<test.transpose();
    for(uint j = 0; j < NO_JOINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < 1e-9);

    //cout<<"\n............................."<<endl;
}

BOOST_AUTO_TEST_CASE(solver_qp_oases_with_equality_constraints)
{
    srand (time(NULL));

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 6;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = false;
    const int NO_WSR = 20;

    // Solve the problem min(||x||), subject Ax=b --> encode the task as constraint
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = I  and g = 0
    // For a 6x6 Constraint matrix this is approx. 3-5 times slower than encoding the task in the cost function as above

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    qp.g.setZero();
    qp.H.setIdentity();
    // Task Jacobian
    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    qp.A = A;
    // Desired task space reference
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;
    qp.b = y;

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    EiquadprogSolver solver;
    solver.setMaxNIter(NO_WSR);

    BOOST_CHECK(solver.getMaxNIter() == NO_WSR);

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);
    //long useconds = end.tv_usec - start.tv_usec;

    /*cout<<"\n----------------------- Test Results ----------------------"<<endl<<endl;
    std::cout<<"Solver took "<<useconds<<" us "<<std::endl;
    cout<<"No of joints: "<<NO_JOINTS<<endl;
    cout<<"No of constraints: "<<NO_CONSTRAINTS<<endl;

    cout<<"\nSolver Input:"<<endl;
    cout<<"Constraint Matrix A:"<<endl; cout<<A<<endl;
    cout<<"Reference: y = "<<y.transpose()<<endl;

    cout<<"\nSolver Output: q_dot = "<<solver_output.transpose()<<endl;*/
    Eigen::VectorXd test = A*solver_output;
    //cout<<"Test: A * q_dot = "<<test.transpose();
    for(uint j = 0; j < y.size(); j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < 1e-9);

    //cout<<"\n............................."<<endl;
}

BOOST_AUTO_TEST_CASE(solver_eiquadprog_bounded)
{
    srand (time(NULL));

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = true;
    const int NO_WSR = 20;

    // Solve the problem min(||Ax-b||) without constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T
    // For a 6x6 Constraint matrix this is approx. 3-5 times faster than encoding the task as constraint as below
    // With warm start, this solver is much faster (approx. 5 times) than in the initial run

    wbc::QuadraticProgram qp;
    qp.resize(NO_JOINTS, NO_EQ_CONSTRAINTS, NO_IN_CONSTRAINTS, WITH_BOUNDS);

    // Task Jacobian
    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    // Desired task space reference
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;

    qp.H = A.transpose()*A;
    qp.g = -(A.transpose()*y).transpose();

    qp.lower_x.setConstant(-0.4);
    qp.upper_x.setConstant(+0.4);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    EiquadprogSolver solver;
    solver.setMaxNIter(NO_WSR);

    BOOST_CHECK(solver.getMaxNIter() == NO_WSR);

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);

    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp.lower_x(j)-1e-9) <= solver_output(j) && solver_output(j) <= (qp.upper_x(j)+1e-9));

}

BOOST_AUTO_TEST_CASE(solver_eiquadprog_unbounded_constraints)
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
    EiquadprogSolver solver_inf, solver_wide;
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
