#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <vector>

enum SimType{FINITE_DIFFERENCE, FINITE_ELEMENT, FINITE_VOLUME};

class Simulation{
public:

	run();

private:
	SimType _sim_type;			// this defines the type of discretization
	Static_Mesh * _mesh;		// this links to the mesh
	std::vector<BoundaryCondition> _boundary_conditions;	// this defines the boudnary conditions
	Equation _equation; 		// this defines the operators and terms involved
	double tstart, tstop;		// this defines the start and stop time

};

#endif _SIMULATION_H