#include "OsqpSolver.hpp"

namespace wbc {

QPSolverRegistry<OsqpSolver> OsqpSolver::reg("osqp");

OsqpSolver::OsqpSolver() : AcadosOcpQPSolver(MODULE_OSQP, "OSQP", 1, 1){
}

}
