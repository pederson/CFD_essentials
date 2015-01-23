#ifndef _CONVERTER_H
#define _CONVERTER_H

#include "Mesh.hpp"
#include "GeometricObject.hpp"
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

// class to define conversion of objects between other objects
// e.g. Hull to PointCloud, or GeometricModel to Mesh, or StaticMesh to MutableMesh
class Converter{
public:

private:

};

// define the functions that build a mesh from a parametric model
Static_Mesh * build_simple_mesh_2d(parametric_model_2d * model,  double res, double xmin, double xmax, double ymin, double ymax, std::vector<double> bg_properties);
//Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);

// helpers
void add_shape_to_mesh(Mutable_Mesh * mesh, geometric_object_2d * shape, double res);
Hull * approximate_parametric_shape_2d(geometric_object_2d * model, double res);


#endif