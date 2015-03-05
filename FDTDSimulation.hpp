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
	void bind_mesh(const RegularMesh & mesh);
	void bind_rel_permittivity(const double * rel_permittivity);
	void bind_rel_permeability(const double * rel_permeability);
	void view_results();
	void output_HDF5(std::string outname="");

	void run(int num_iters = -1);

private:

	void preRunCheck();
	void allocate_fields();
	void allocate_PML();
	void allocate_simdata();

	const RegularMesh * _mesh;

	// internal variables
	const double _eps0 = 8.854e-12;
	const double _mu0 = 4*3.14159e-7;
	const double _c0 = 2.99792458e+8;
	double _dx, _dt, _CourantFactor;
	double _nPML, _num_iters, _current_iter, _tcur;

	// user options
	bool _output_HDF5;
	std::string _output_HDF5_name;
	bool _visualize_results;

	const double * _rel_permittivity;
	const double * _rel_permeability;

	// field stuff
	double * D_z;
	double * I_Hx;
	double * I_Hy;
	double * E_z;
	double * H_y;
	double * H_x;
	double * gaz;

	// PML stuff
	double * gi2;
	double * gi3;
	double * gj2;
	double * gj3;
	double * fi1;
	double * fi2;
	double * fi3;
	double * fj1;
	double * fj2;
	double * fj3;



};

#endif