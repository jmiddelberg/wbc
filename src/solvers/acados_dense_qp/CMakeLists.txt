SET(TARGET_NAME wbc-acados_dense_qp)

# Common base of the solvers that wbc calls through the dense QP interface of acados (DAQP, ProxQP,
# qpOASES). acados is found in ../CMakeLists.txt.

set(SOURCES AcadosDenseQPSolver.cpp)
set(HEADERS AcadosDenseQPSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-core)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
# acados is private: the header has no acados types, so code that uses the solvers needs neither
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
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/acados_dense_qp)
