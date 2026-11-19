SET(TARGET_NAME wbc-solvers-clarabel)

# Clarabel comes with acados (ACADOS_WITH_CLARABEL) and is called through its OCP QP solver module,
# see ../acados_qp. acados is found in ../CMakeLists.txt.
wbc_require_acados_solver(ACADOS_WITH_CLARABEL)

set(SOURCES ClarabelSolver.cpp)
set(HEADERS ClarabelSolver.hpp)

list(APPEND PKGCONFIG_REQUIRES wbc-acados_qp)
string (REPLACE ";" " " PKGCONFIG_REQUIRES "${PKGCONFIG_REQUIRES}")

add_library(${TARGET_NAME} SHARED ${SOURCES} ${HEADERS})
target_link_libraries(${TARGET_NAME}
                      PUBLIC  wbc-acados_qp)

set_target_properties(${TARGET_NAME} PROPERTIES
       VERSION ${PROJECT_VERSION}
       SOVERSION ${API_VERSION})

install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION lib)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}.pc DESTINATION lib/pkgconfig)
INSTALL(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME}/solvers/clarabel)

add_subdirectory(test)
