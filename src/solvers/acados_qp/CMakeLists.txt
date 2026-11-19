SET(TARGET_NAME wbc-acados_qp)

# Common base of the solvers that wbc calls through acados: the dense QP interface (DAQP, ProxQP,
# qpOASES) and the OCP QP solver modules (OSQP, Clarabel). acados is found in ../CMakeLists.txt.

set(SOURCES AcadosQPSolver.cpp
            AcadosDenseQPSolver.cpp
            AcadosOcpQPSolver.cpp)
set(HEADERS AcadosQPSolver.hpp
            AcadosDenseQPSolver.hpp
            AcadosOcpQPSolver.hpp)
# AcadosQPBackend.hpp is internal to this library and not installed

list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS} AcadosQPBackend.hpp)
# acados is private: the headers have no acados types, so code that uses the solvers needs neither
# the acados include paths nor links acados itself
target_link_libraries(${TARGET_NAME}
                      PUBLIC  wbc-core
                      PRIVATE acados)

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/acados_qp)
