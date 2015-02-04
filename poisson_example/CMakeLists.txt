CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(visualixer)

SET(CMAKE_CXX_FLAGS "-std=c++11 -w -Wall -Werror")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

FIND_PACKAGE( PkgConfig REQUIRED )
PKG_SEARCH_MODULE( GLFW3 REQUIRED glfw3 )

FIND_PACKAGE(PETSc REQUIRED)

INCLUDE_DIRECTORIES( ${GLFW3_INCLUDE_DIR} )
ADD_EXECUTABLE( poisson_test PoissonSolve.cpp ../Visualixer.cpp ../VisualixerMesh.cpp ../PointCloud.cpp ../Mesh.cpp ../GeometricObject.cpp ../ColorRamp.cpp ../Converter.cpp ../Hull.cpp)

FOREACH(item ${GLFW3_STATIC_LIBRARIES})
    MESSAGE(STATUS "  using lib: " ${item})
ENDFOREACH()

TARGET_LINK_LIBRARIES( poisson_test GLEW glfw3 ${GLFW3_STATIC_LIBRARIES} ${OPENGL_gl_LIBRARY} ${PETSC_LIBRARIES})
