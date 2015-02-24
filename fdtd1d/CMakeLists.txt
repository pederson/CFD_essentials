# you can compile this with the following command:
# cmake .. -DPETSC_DIR=/home/dmp942/GitProjects/petsc/ -DPETSC_ARCH=arch-linux2-c-debug

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(fdtd1d)

SET(CMAKE_CXX_COMPILER "${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx")
#SET(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
#SET(CMAKE_CXX_FLAGS "-std=c++11 -w -Wall -Werror -fopenmp")
SET(CMAKE_CXX_FLAGS "-std=c++11 -w -Wall -Werror -fopenmp")
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# PETSc
#SET(PETSC_DIR "/home/dmp942/GitProjects/petsc")
#SET(PETSC_ARCH "arch-linux2-c-debug")
FIND_PACKAGE(PETSc REQUIRED)
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

# package config
FIND_PACKAGE( PkgConfig REQUIRED )

# GLFW
PKG_SEARCH_MODULE( GLFW3 REQUIRED glfw3 )
INCLUDE_DIRECTORIES(${GLFW3_INCLUDE_DIR} )

# HDF5
FIND_PACKAGE(HDF5 NAMES hdf5)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})



ADD_EXECUTABLE( fdtd1d fdtd1d.cpp ../Visualixer.cpp ../VisualixerMesh.cpp ../Mesh.cpp ../RegularMesh.cpp ../ColorRamp.cpp ../SimulationData.cpp)

FOREACH(item ${GLFW3_STATIC_LIBRARIES})
    MESSAGE(STATUS "  using lib: " ${item})
ENDFOREACH()

TARGET_LINK_LIBRARIES( fdtd1d GLEW glfw3 ${GLFW3_STATIC_LIBRARIES} ${OPENGL_gl_LIBRARY} ${PETSC_LIBRARIES} ${HDF5_LIBRARIES})
