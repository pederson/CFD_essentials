#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <vector>

enum SimType{FINITE_DIFFERENCE, FINITE_ELEMENT, FINITE_VOLUME};

class Simulation{
public:

	run();

private:
	SimType _sim_type;
	Static_Mesh * _mesh;
	std::vector<BoundaryCondition> _boundary_conditions;
	Equation _equation;

};

#endif _SIMULATION_H