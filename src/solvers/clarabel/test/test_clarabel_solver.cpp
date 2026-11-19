#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "solvers/clarabel/ClarabelSolver.hpp"
#include "solvers/acados_qp/test/QPTestCases.hpp"

using namespace wbc;

// Clarabel is an interior-point method: the solution is as accurate as the tolerances that setTolerance() sets,
// up to the conditioning of the problem
static const double TOL = 1e-6;

static std::unique_ptr<AcadosQPSolver> makeSolver(){
    std::unique_ptr<AcadosQPSolver> solver(new ClarabelSolver());
    solver->setTolerance(1e-9);
    return solver;
}

WBC_ACADOS_QP_TEST_CASES(solver_clarabel, makeSolver, TOL)

BOOST_AUTO_TEST_CASE(solver_clarabel_options)
{
    // Clarabel does not warm start
    acados_qp_test::commonOptions([](){ return std::unique_ptr<AcadosQPSolver>(new ClarabelSolver()); }, 0, 0);
}

BOOST_AUTO_TEST_CASE(solver_clarabel_plugin_name)
{
    std::unique_ptr<QPSolver> solver(QPSolverFactory::createInstance("clarabel"));
    BOOST_CHECK(dynamic_cast<ClarabelSolver*>(solver.get()) != nullptr);
}
