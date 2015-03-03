#ifndef _FDTDSIMULATION_H
#define _FDTDSIMULATION_H

#include "Simulation.hpp"

#include <vector>
#include <string>

enum BoundaryLocation{BOUNDARY_RIGHT, BOUNDARY_LEFT, BOUNDARY_TOP, BOUNDARY_BOTTOM, BOUNDARY_FRONT, BOUNDARY_BACK};
enum BoundaryCondition{BOUNDARY_PML, BOUNDARY_PERIODIC, BOUNDARY_PEC, BOUNDARY_PMC};

class FDTDSimulation : public Simulation{
public:

	// inspectors

	// mutators
	void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void add_gaussian_source(double t0, double spread, double xloc, double yloc=0, double zloc=0);
	void add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc=0, double zloc=0);
	void bind_rel_permittivity(const double * rel_permittivity);
	void set_view_results(bool opt);
	void set_output_HDF5(bool opt, std::string outname);

	void run();

private:

	void preRunCheck();

	bool _output_HDF5;
	std::string _output_HDF5_name;
	bool _visualize_results;

	const double * _rel_permittivity;



};

#endif