#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "solvers/osqp/OsqpSolver.hpp"
#include "solvers/acados_qp/test/QPTestCases.hpp"

using namespace wbc;

// OSQP stops when the residuals are below eps_abs + eps_rel * (size of the terms), which setTolerance() sets. The
// polished solution is usually more accurate.
static const double TOL = 1e-6;

static std::unique_ptr<AcadosQPSolver> makeSolver(){
    std::unique_ptr<AcadosQPSolver> solver(new OsqpSolver());
    solver->setTolerance(1e-8);
    return solver;
}

WBC_ACADOS_QP_TEST_CASES(solver_osqp, makeSolver, TOL)

BOOST_AUTO_TEST_CASE(solver_osqp_options)
{
    acados_qp_test::commonOptions([](){ return std::unique_ptr<AcadosQPSolver>(new OsqpSolver()); }, 1, 1);
}

BOOST_AUTO_TEST_CASE(solver_osqp_plugin_name)
{
    std::unique_ptr<QPSolver> solver(QPSolverFactory::createInstance("osqp"));
    BOOST_CHECK(dynamic_cast<OsqpSolver*>(solver.get()) != nullptr);
}
