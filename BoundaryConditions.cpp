#include "BoundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(Static_Mesh * mesh, BC_Type type, std::string property, double val){
	_mesh = mesh;
	_type = type;
	_associated_property = property;
	_associated_value = val;
}

BoundaryConditions::BoundaryConditions(){
	_mesh = NULL;
}

BoundaryConditions::~BoundaryConditions(){

}
