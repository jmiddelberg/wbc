#include "ClarabelSolver.hpp"

namespace wbc {

QPSolverRegistry<ClarabelSolver> ClarabelSolver::reg("clarabel");

ClarabelSolver::ClarabelSolver() : AcadosOcpQPSolver(MODULE_CLARABEL, "Clarabel", 0, 0){
}

}
