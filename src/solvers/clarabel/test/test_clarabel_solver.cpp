#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <vector>
#include "../../../core/QuadraticProgram.hpp"
#include "../ClarabelSolver.hpp"

using namespace wbc;
using namespace std;

// Tolerance for checking the solution. Clarabel is an interior-point solver, so its default
// termination tolerances (~1e-8) are looser than the KKT-exact solvers.
static const double TOL = 1e-5;

BOOST_AUTO_TEST_CASE(solver_clarabel_without_constraints)
{
    srand (time(NULL));

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = false;

    // Solve the problem min(||Ax-b||) without constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T

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

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    ClarabelSolver solver;

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);

    Eigen::VectorXd test = A*solver_output;
    //cout<<"Test: A * q_dot = "<<test.transpose();
    for(uint j = 0; j < NO_JOINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);

    //cout<<"\n............................."<<endl;
}

BOOST_AUTO_TEST_CASE(solver_clarabel_with_equality_constraints)
{
    srand (time(NULL));

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

    ClarabelSolver solver;

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);
    //long useconds = end.tv_usec - start.tv_usec;

    Eigen::VectorXd test = A*solver_output;
    //cout<<"Test: A * q_dot = "<<test.transpose();
    for(uint j = 0; j < NO_EQ_CONSTRAINTS; j++)
        BOOST_CHECK(fabs(test(j) - y(j)) < TOL);

    //cout<<"\n............................."<<endl;
}

BOOST_AUTO_TEST_CASE(solver_clarabel_with_inequality_constraints)
{
    srand (time(NULL));

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
    // Task Jacobian
    Eigen::MatrixXd A(6,6);
    A << 0.642, 0.706, 0.565,  0.48,  0.59, 0.917,
         0.553, 0.087,  0.43,  0.71, 0.148,  0.87,
         0.249, 0.632, 0.711,  0.13, 0.426, 0.963,
         0.682, 0.123, 0.998, 0.716, 0.961, 0.901,
         0.891, 0.019, 0.716, 0.534, 0.725, 0.633,
         0.315, 0.551, 0.462, 0.221, 0.638, 0.244;
    qp.C = A;
    // Desired task space reference
    Eigen::VectorXd y(6);
    y << 0.833, 0.096, 0.078, 0.971, 0.883, 0.366;
    qp.lower_y = y - Eigen::VectorXd::Constant(qp.nq, 1e-1);
    qp.upper_y = y + Eigen::VectorXd::Constant(qp.nq, 1e-1);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    ClarabelSolver solver;

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);
    //long useconds = end.tv_usec - start.tv_usec;

    Eigen::VectorXd test = A*solver_output;

    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp.lower_y(j)-TOL) <= test(j) && test(j) <= (qp.upper_y(j)+TOL));
}

BOOST_AUTO_TEST_CASE(solver_clarabel_bounded)
{
    srand (time(NULL));

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 0;
    const bool WITH_BOUNDS = true;

    // Solve the problem min(||Ax-b||) with bound constraints --> encode the task as part of the cost function
    // Standard form of QP is x^T*H*x + x^T*g --> Choose H = A^T*A and g = -(A^T*y)^T

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

    qp.lower_x.setConstant(-1e10);
    qp.upper_x.setConstant(+1e10);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    ClarabelSolver solver;

    Eigen::VectorXd solver_output;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    BOOST_CHECK_NO_THROW(solver.solve(hqp, solver_output));
    gettimeofday(&end, NULL);

    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp.lower_x(j)-TOL) <= solver_output(j) && solver_output(j) <= (qp.upper_x(j)+TOL));

}

BOOST_AUTO_TEST_CASE(solver_clarabel_unbounded_constraints)
{
    // Bounds marked with wbc::INF are absent. ClarabelSolver replaces their right-hand side with an
    // actual infinity, which makes Clarabel's presolve delete the row. That must not change the
    // solution: solving the same problem with bounds that are finite but far too wide to ever
    // become active has to give the same answer, in no more iterations.

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 2;
    const int NO_IN_CONSTRAINTS = 6;

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

    // The 6 unbounded inequality rows and the 6 rows of the three free variables are dropped
    // before Clarabel sees the problem, the 6 rows of the bounded variables and the ZeroCone of
    // the equalities are kept.
    ClarabelSolver solver_inf, solver_wide;
    BOOST_CHECK_EQUAL(solver_inf.getInfinity(), wbc::INF);
    Eigen::VectorXd out_inf, out_wide;
    BOOST_CHECK_NO_THROW(solver_inf.solve(hqp_inf, out_inf));
    BOOST_CHECK_NO_THROW(solver_wide.solve(hqp_wide, out_wide));

    // Same problem, so the same solution ...
    BOOST_CHECK_EQUAL(out_inf.size(), NO_JOINTS);
    for(uint j = 0; j < NO_JOINTS; ++j)
        BOOST_CHECK_SMALL(out_inf(j) - out_wide(j), TOL);

    // ... reached in no more iterations than the formulation that keeps all rows needs.
    BOOST_CHECK_LE(solver_inf.getNter(), solver_wide.getNter());

    // The bounds that are not marked as absent must still be enforced ...
    for(uint j = 3; j < NO_JOINTS; ++j)
        BOOST_CHECK((qp_inf.lower_x(j)-TOL) <= out_inf(j) && out_inf(j) <= (qp_inf.upper_x(j)+TOL));
    // ... and at least one of them has to be active, otherwise dropping every bound would pass too
    BOOST_CHECK(out_inf.tail(3).cwiseAbs().maxCoeff() > 0.4 - 1e-2);

    // The equality rows are never touched by the substitution
    BOOST_CHECK_SMALL((qp_inf.A*out_inf - qp_inf.b).cwiseAbs().maxCoeff(), TOL);
}

BOOST_AUTO_TEST_CASE(solver_clarabel_linear_solver_selection)
{
    // The linear system solver that Clarabel uses to factorize the KKT system can be selected by
    // enum or by name, see https://clarabel.org/stable/user_guide_linsolvers/. Which methods are
    // available depends on the cargo features libclarabel_c was built with, so this test runs over
    // whatever this build offers and requires that all of them solve the same problem alike.

    const int NO_JOINTS = 6;
    const int NO_EQ_CONSTRAINTS = 0;
    const int NO_IN_CONSTRAINTS = 6;
    const bool WITH_BOUNDS = true;

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
    // One-sided inequalities and bounds that the unconstrained solution violates, so that the
    // problem is feasible (x = 0 is) but all cones are actually in play
    qp.C = A;
    qp.lower_y.setConstant(-wbc::INF);
    qp.upper_y = y - Eigen::VectorXd::Constant(NO_IN_CONSTRAINTS, 5e-2);
    qp.lower_x.setConstant(-1.0);
    qp.upper_x.setConstant(+1.0);

    BOOST_CHECK(qp.isValid());
    wbc::HierarchicalQP hqp;
    hqp << qp;

    // auto and qdldl are always compiled in, the names round trip
    const std::vector<std::string> available = ClarabelSolver::availableLinearSolvers();
    BOOST_CHECK(std::find(available.begin(), available.end(), "auto") != available.end());
    BOOST_CHECK(std::find(available.begin(), available.end(), "qdldl") != available.end());
    for(const std::string& name : available)
        BOOST_CHECK_EQUAL(ClarabelSolver::linearSolverName(ClarabelSolver::linearSolverFromName(name)), name);
    cout<<"Available Clarabel linear system solvers:";
    for(const std::string& name : available)
        cout<<" "<<name;
    cout<<endl;

    // wbc defaults to qdldl, and auto has to resolve to one of the concrete methods
    ClarabelSolver solver;
    BOOST_CHECK_EQUAL(solver.getLinearSolver(), clarabel_linsolver_qdldl);
    Eigen::VectorXd reference;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, reference));
    BOOST_CHECK_EQUAL(solver.getLinearSolverUsed(), clarabel_linsolver_qdldl);

    ClarabelSolver solver_auto;
    solver_auto.setLinearSolver(clarabel_linsolver_auto);
    Eigen::VectorXd out_auto;
    BOOST_CHECK_NO_THROW(solver_auto.solve(hqp, out_auto));
    BOOST_CHECK(solver_auto.getLinearSolverUsed() != clarabel_linsolver_auto);
    BOOST_CHECK(ClarabelSolver::linearSolverAvailable(solver_auto.getLinearSolverUsed()));

    // Every available method has to be selectable, has to be the one Clarabel then actually uses,
    // and has to end up at the same solution
    for(const std::string& name : available){
        const ClarabelLinearSolver method = ClarabelSolver::linearSolverFromName(name);

        ClarabelSolver solver_ls;
        BOOST_CHECK_NO_THROW(solver_ls.setLinearSolver(name));
        BOOST_CHECK_EQUAL(solver_ls.getLinearSolver(), method);
        // ClarabelLinearSolver is numbered like the enum Clarabel expects over the FFI
        BOOST_CHECK_EQUAL(static_cast<int>(solver_ls.getOptions().direct_solve_method), static_cast<int>(method));

        Eigen::VectorXd out;
        BOOST_CHECK_NO_THROW(solver_ls.solve(hqp, out));
        if(method != clarabel_linsolver_auto)
            BOOST_CHECK_EQUAL(solver_ls.getLinearSolverUsed(), method);

        for(uint j = 0; j < NO_JOINTS; ++j)
            BOOST_CHECK_SMALL(out(j) - reference(j), TOL);

        // Setting the same method by enum has to be equivalent to setting it by name
        ClarabelSolver solver_enum;
        BOOST_CHECK_NO_THROW(solver_enum.setLinearSolver(method));
        BOOST_CHECK_EQUAL(solver_enum.getLinearSolver(), method);
    }

    // Methods that are not part of this build are rejected right away, and so are unknown names
    const ClarabelLinearSolver all[] = {clarabel_linsolver_auto,
                                        clarabel_linsolver_qdldl,
                                        clarabel_linsolver_faer,
                                        clarabel_linsolver_pardiso_mkl,
                                        clarabel_linsolver_pardiso_panua};
    for(ClarabelLinearSolver method : all){
        if(!ClarabelSolver::linearSolverAvailable(method)){
            BOOST_CHECK_THROW(solver.setLinearSolver(method), std::runtime_error);
            BOOST_CHECK_THROW(solver.setLinearSolver(ClarabelSolver::linearSolverName(method)), std::runtime_error);
        }
    }
    BOOST_CHECK_THROW(ClarabelSolver::linearSolverFromName("cholmod"), std::invalid_argument);
    BOOST_CHECK_THROW(solver.setLinearSolver(std::string("cholmod")), std::invalid_argument);

    // A rejected method must not have changed anything
    BOOST_CHECK_EQUAL(solver.getLinearSolver(), clarabel_linsolver_qdldl);

    // setOptions() replaces the whole settings struct, including the linear system solver
    clarabel::DefaultSettings<double> opt = clarabel::DefaultSettings<double>::default_settings();
    opt.verbose = false;
    opt.direct_solve_method = clarabel::ClarabelDirectSolveMethods::QDLDL;
    solver.setOptions(opt);
    BOOST_CHECK_EQUAL(solver.getLinearSolver(), clarabel_linsolver_qdldl);
    Eigen::VectorXd out_qdldl;
    BOOST_CHECK_NO_THROW(solver.solve(hqp, out_qdldl));
    BOOST_CHECK_EQUAL(solver.getLinearSolverUsed(), clarabel_linsolver_qdldl);
}
