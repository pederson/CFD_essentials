#ifndef _BOUNDARYCONDITIONS_H
#define _BOUNDARYCONDITIONS_H

#include "Mesh.hpp"

#include <vector>
#include <string>
#include <iostream>

// this class should allow for definition of regular grid and unstructured grid 
// boundary conditions. 
// This could be Dirichlet, Neumann, Cauchy, periodic, or other possibilities. 
// They could also be static or dynamic
// They could be varying or constant along the boundary
enum BC_Type{BC_DIRICHLET, BC_NEUMANN, BC_CAUCHY, PERIODIC}
class BoundaryConditions{
public:
	BoundaryConditions(Static_Mesh * mesh, BC_Type type, std::string property, double val);
	BoundaryConditions();
	~BoundaryConditions();

	void set_associated_mesh(Static_Mesh * mesh){_mesh = mesh;};
	void set_BC_type(BC_Type type){_type = type;};
	void set_boundary_property(std::string prop_name){_associated_property = prop_name;};
	void set_boundary_value(double val){_associated_value = val;};
	void set_boundary_inds(std::vector<unsigned int> inds){_boundary_inds = inds;};

private:
	Static_Mesh * _mesh;
	BC_Type _type;
	std::string _associated_property;
	double _associated_value;

	std::vector<unsigned int> _boundary_inds;
};


#endif