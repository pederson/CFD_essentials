#ifndef _MODEL2MESH_
#define _MODEL2MESH_

#include "mesh_class.hpp"
#include "geometric_object.hpp"

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

// define the functions that build a mesh from a parametric model
Mesh * build_simple_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);
//Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);




#endif