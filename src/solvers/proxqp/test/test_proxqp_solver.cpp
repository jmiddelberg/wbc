#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "solvers/proxqp/ProxQPSolver.hpp"
#include "solvers/acados_dense_qp/test/DenseQPTestCases.hpp"

using namespace wbc;

// ProxQP stops when the primal and dual residuals are below eps_abs, which setTolerance() sets
static const double TOL = 1e-8;

static std::unique_ptr<AcadosDenseQPSolver> makeSolver(){
    std::unique_ptr<AcadosDenseQPSolver> solver(new ProxQPSolver());
    solver->setTolerance(1e-9);
    return solver;
}

WBC_DENSE_QP_TEST_CASES(solver_proxqp, makeSolver, TOL)

BOOST_AUTO_TEST_CASE(solver_proxqp_options)
{
    dense_qp_test::commonOptions([](){ return std::unique_ptr<AcadosDenseQPSolver>(new ProxQPSolver()); }, 2, 3);
}

BOOST_AUTO_TEST_CASE(solver_proxqp_dense_backends)
{
    // Every backend has to give the same solution. The constraints are built around a feasible
    // point, the cost pulls the solution away from it, so some constraints are active.
    wbc::QuadraticProgram qp;
    qp.resize(6, 2, 6, true);
    Eigen::MatrixXd A = dense_qp_test::taskJacobian();
    Eigen::VectorXd y = dense_qp_test::taskReference();
    const Eigen::VectorXd x_feasible = Eigen::VectorXd::Constant(6, 0.1);
    qp.H = A.transpose()*A + 1e-6*Eigen::MatrixXd::Identity(6,6);
    qp.g = -(A.transpose()*(2*y));
    qp.A = A.topRows(2);
    qp.b = qp.A*x_feasible;
    qp.C = A;
    qp.lower_y.setConstant(-wbc::INF);
    qp.upper_y = A*x_feasible + Eigen::VectorXd::Constant(6, 0.05);
    qp.lower_x.setConstant(-0.4);
    qp.upper_x.setConstant(+0.4);

    ProxQPSolver solver;
    solver.setTolerance(1e-9);
    BOOST_CHECK_EQUAL(solver.getDenseBackend(), "AUTOMATIC");
    Eigen::VectorXd x_automatic;
    dense_qp_test::solve(solver, qp, x_automatic);

    for(const std::string backend : {"PRIMAL_DUAL_LDLT", "PRIMAL_LDLT", "AUTOMATIC"}){
        // a change after the first solve takes effect at the next solve
        BOOST_CHECK_NO_THROW(solver.setDenseBackend(backend));
        Eigen::VectorXd x;
        if(dense_qp_test::solve(solver, qp, x))
            BOOST_CHECK_MESSAGE((x - x_automatic).cwiseAbs().maxCoeff() < 1e-6, backend << ": solution differs by " << (x - x_automatic).cwiseAbs().maxCoeff());
    }

    BOOST_CHECK_THROW(solver.setDenseBackend("SPARSE"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(solver_proxqp_plugin_name)
{
    std::unique_ptr<QPSolver> solver(QPSolverFactory::createInstance("proxqp"));
    BOOST_CHECK(dynamic_cast<ProxQPSolver*>(solver.get()) != nullptr);
}
