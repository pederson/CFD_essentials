CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(visualixer)

SET(CMAKE_CXX_FLAGS "-std=c++11 -w -Wall -Werror")

FIND_PACKAGE( PkgConfig REQUIRED )
PKG_SEARCH_MODULE( GLFW3 REQUIRED glfw3 )

INCLUDE_DIRECTORIES( ${GLFW3_INCLUDE_DIR} )
ADD_EXECUTABLE( visualixer_test visualixer.cpp PointCloud.cpp mesh_class.cpp geometric_object.cpp colorramp.cpp)

FOREACH(item ${GLFW3_STATIC_LIBRARIES})
    MESSAGE(STATUS "  using lib: " ${item})
ENDFOREACH()

TARGET_LINK_LIBRARIES( visualixer_test GLEW glfw3 ${GLFW3_STATIC_LIBRARIES} ${OPENGL_gl_LIBRARY})
