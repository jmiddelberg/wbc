#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "solvers/qpoases/QPOasesSolver.hpp"
#include "solvers/acados_qp/test/QPTestCases.hpp"

using namespace wbc;

// qpOASES is an active-set method, so the solution is exact up to the conditioning of the problem
static const double TOL = 1e-9;

static std::unique_ptr<AcadosDenseQPSolver> makeSolver(){
    return std::unique_ptr<AcadosDenseQPSolver>(new QPOASESSolver());
}

WBC_ACADOS_QP_TEST_CASES(solver_qpoases, makeSolver, TOL)

BOOST_AUTO_TEST_CASE(solver_qpoases_options)
{
    // warm start 2 is a hot start, which is only valid for constant matrices
    acados_qp_test::commonOptions(makeSolver, 1, 1);
}

BOOST_AUTO_TEST_CASE(solver_qpoases_plugin_name)
{
    std::unique_ptr<QPSolver> solver(QPSolverFactory::createInstance("qpoases"));
    BOOST_CHECK(dynamic_cast<QPOASESSolver*>(solver.get()) != nullptr);
}
