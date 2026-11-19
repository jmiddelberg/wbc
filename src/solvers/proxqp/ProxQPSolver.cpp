#include "ProxQPSolver.hpp"

#include <acados_c/dense_qp_interface.h>
#include <stdexcept>

namespace wbc {

QPSolverRegistry<ProxQPSolver> ProxQPSolver::reg("proxqp");

// Warm start level 2 keeps the proximal parameters, as this wrapper did when it called ProxQP directly
// (WARM_START_WITH_PREVIOUS_RESULT). acados accepts levels up to 3.
ProxQPSolver::ProxQPSolver() : AcadosDenseQPSolver(DENSE_QP_PROXQP, "ProxQP", 2, 3){
}

/// Validated here, because acados terminates the process on a backend it does not know
void ProxQPSolver::setDenseBackend(const std::string& backend){
    if(backend != "AUTOMATIC" && backend != "PRIMAL_DUAL_LDLT" && backend != "PRIMAL_LDLT")
        throw std::invalid_argument("Invalid ProxQP dense backend '" + backend + "', must be one of AUTOMATIC, PRIMAL_DUAL_LDLT, PRIMAL_LDLT");
    dense_backend = backend;
    reset();
}

void ProxQPSolver::setSolverOptions(){
    setOption("proxqp_dense_backend", (void*)dense_backend.c_str());
}

} // namespace wbc
