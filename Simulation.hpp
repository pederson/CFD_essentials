#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "BoundaryConditions.hpp"
#include "InitialConditions.hpp"
#include "Equation.hpp"
#include "Mesh.hpp"
#include "LinAlgWrapper.hpp"

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
	void bind_mesh(const Mesh & mesh) {_mesh = &mesh;};
	//void add_equation(const Equation & eq) {_equations.push_back(eq);};
	void add_boundary_condition(const BoundaryCondition & BC) {_boundary_conditions.push_back(BC);};
	void add_initial_condition(const InitialCondition & IC) {_initial_conditions.push_back(IC);};
	void set_time_span(double tstart, double dt, double tstop) {_tstart=tstart; _tstop=tstop; _dt=dt;};
	//set_solution_method

	virtual void run();

private:

	virtual void preRunCheck();

	SimType _sim_type;			// this defines the type of discretization
	const Mesh * _basemesh;		// this links to the mesh
	SimulationData _simdata;	// this holds the data that comes from the mesh
	std::vector<BoundaryCondition> _boundary_conditions;	// this defines the boudnary conditions
	std::vector<InitialCondition> _initial_conditions;
	std::vector<Equation> _equations; 		// this defines the equations involved
	double _tstart, _dt, _tstop;		// this defines the start and stop time

	// define linear algebra solution method
	//LinVector _vec;


};


#endif