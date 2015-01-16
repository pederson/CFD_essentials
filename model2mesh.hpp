#ifndef _MODEL2MESH_H
#define _MODEL2MESH_H

#include "mesh_class.hpp"
#include "geometric_object.hpp"
#include "Hull.hpp"

#include <iostream>
#include <vector>
#include <string>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

// define the functions that build a mesh from a parametric model
Mesh * build_simple_mesh_2d(parametric_model_2d * model,  double res, double xmin, double xmax, double ymin, double ymax, std::vector<double> bg_properties);
//Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);

// helpers
void add_shape_to_mesh(Mesh * meshmodel, geometric_object_2d * shape, double res);
Hull * approximate_parametric_shape_2d(geometric_object_2d * model, double res);


#endif