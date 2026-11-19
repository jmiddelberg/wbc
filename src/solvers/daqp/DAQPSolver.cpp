#include "DAQPSolver.hpp"

#include <acados_c/dense_qp_interface.h>

namespace wbc {

QPSolverRegistry<DAQPSolver> DAQPSolver::reg("daqp");

DAQPSolver::DAQPSolver() : AcadosDenseQPSolver(DENSE_QP_DAQP, "DAQP", 1, 1){
}

} // namespace wbc
