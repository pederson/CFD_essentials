#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "BoundaryConditions.hpp"
#include "InitialConditions.hpp"
#include "Equation.hpp"
#include "Mesh.hpp"

#include <vector>

enum SimType{SIM_FD, SIM_FDTD, SIM_FV, SIM_FVTD, SIM_FE, SIM_FETD};

// this class should gather together equations, boundary conditions, and a mesh
// and perform checking on them before running a simulation
class Simulation{
public:
	Simulation(SimType type, const Mesh & mesh);
	~Simulation();

	// inspectors
	void print_summary() const;

	// mutators
	bind_mesh(const Mesh & mesh);
	add_equation(const Equation & eq) {_equations.push_back(eq);};
	add_boundary_condition(const BoundaryCondition & BC) {_boundary_conditions.push_back(BC);};
	add_initial_condition(const InitialCondition & IC) {_initial_conditions.push_back(IC);};

	run();

private:
	SimType _sim_type;			// this defines the type of discretization
	const Mesh * _mesh;		// this links to the mesh
	SimulationData _simdata;	// this holds the data that comes from the mesh
	std::vector<BoundaryCondition> _boundary_conditions;	// this defines the boudnary conditions
	std::vector<InitialCondition> _initial_conditions;
	std::vector<Equation> _equations; 		// this defines the equations involved
	double tstart, tstop;		// this defines the start and stop time

	// define linear algebra solution method
	LinVector


};

class FDTDSimulation : public Simulation{
public:

private:

};

#endif _SIMULATION_H