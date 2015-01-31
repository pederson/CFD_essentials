#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "BoundaryConditions.hpp"
#include "Mesh.hpp"
#include ""

#include <vector>

enum SimType{FINITE_DIFFERENCE, FINITE_ELEMENT, FINITE_VOLUME};

// this class should gather together equations, boundary conditions, and a mesh
// and perform checking on them before running a simulation
class Simulation{
public:
	Simulation(SimType type, Static_Mesh * mesh);
	~Simulation();

	bind_mesh(Static_Mesh * mesh);
	add_equation(Equation eq);

	run();

private:
	SimType _sim_type;			// this defines the type of discretization
	Static_Mesh * _mesh;		// this links to the mesh
	SimulationData _simdata;	// this holds the data that comes from the mesh
	std::vector<BoundaryCondition> _boundary_conditions;	// this defines the boudnary conditions
	std::vector<Equation> _equations; 		// this defines the equations involved
	double tstart, tstop;		// this defines the start and stop time

	// define linear algebra solution method


};

#endif _SIMULATION_H