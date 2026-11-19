#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "solvers/daqp/DAQPSolver.hpp"
#include "solvers/acados_qp/test/QPTestCases.hpp"

using namespace wbc;

// DAQP terminates in a KKT point, so the solution is exact up to the conditioning of the problem
static const double TOL = 1e-9;

static std::unique_ptr<AcadosDenseQPSolver> makeSolver(){
    return std::unique_ptr<AcadosDenseQPSolver>(new DAQPSolver());
}

WBC_ACADOS_QP_TEST_CASES(solver_daqp, makeSolver, TOL)

BOOST_AUTO_TEST_CASE(solver_daqp_options)
{
    // warm start 2 keeps the factorization of the previous call, which is wrong for changing matrices
    acados_qp_test::commonOptions(makeSolver, 1, 1);
}

BOOST_AUTO_TEST_CASE(solver_daqp_plugin_name)
{
    std::unique_ptr<QPSolver> solver(QPSolverFactory::createInstance("daqp"));
    BOOST_CHECK(dynamic_cast<DAQPSolver*>(solver.get()) != nullptr);
}
