#include "QPOasesSolver.hpp"

#include <acados_c/dense_qp_interface.h>

namespace wbc{

QPSolverRegistry<QPOASESSolver> QPOASESSolver::reg("qpoases");

QPOASESSolver::QPOASESSolver() : AcadosDenseQPSolver(DENSE_QP_QPOASES, "qpOASES", 1, 1){
}

}
