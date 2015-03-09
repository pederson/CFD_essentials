#ifndef _FDTDSIMULATION_H
#define _FDTDSIMULATION_H

//#include "Simulation.hpp"
#include "RegularMesh.hpp"
#include "VisualixerSimulation.hpp"
#include "SignalGenerator.hpp"

#include <vector>
#include <string>
#include <iostream>

enum BoundaryLocation{BOUNDARY_RIGHT, BOUNDARY_LEFT, BOUNDARY_TOP, BOUNDARY_BOTTOM, BOUNDARY_FRONT, BOUNDARY_BACK};
enum BoundaryCondition{BOUNDARY_PML, BOUNDARY_PERIODIC, BOUNDARY_PEC, BOUNDARY_PMC};

class FDTDSimulation{
public:

	FDTDSimulation();
	~FDTDSimulation();

	// inspectors

	// mutators
	void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void add_gaussian_source(double t0, double spread, double xloc, double yloc=0, double zloc=0);
	void add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc=0, double zloc=0);
	void add_gaussian_modulator(unsigned int signal_idx, double t0, double spread);
	void add_sinusoidal_modulator(unsigned int signal_idx, double freq_Hz, double phase=0.0);
	void bind_mesh(const RegularMesh & mesh);
	void bind_rel_permittivity(const double * rel_permittivity);
	//void bind_rel_permeability(const double * rel_permeability);
	void bind_current_density_x(const double * current_density_x);
	void bind_current_density_x(const double * velocity_x, const double * density_multiplier_x, const double cdxmult=1.0);
	void bind_current_density_y(const double * current_density_y);
	void bind_current_density_y(const double * velocity_y, const double * density_multiplier_y, const double cdymult=1.0);
	void bind_current_density_z(const double * current_density_z);
	void bind_current_density_z(const double * velocity_z, const double * density_multiplier_z, const double cdzmult=1.0);
	void set_num_iters(unsigned int num_iters) {_num_iters = num_iters;};
	//void set_invariant_time(unsigned int timespan);
	//void set_invariant_resolution(unsigned int res);
	void view_results();
	void output_HDF5(std::string outname="");

	void run(int num_iters = -1);

private:

	void preRunCheck();
	void allocate_fields();
	void allocate_PML();
	void allocate_simdata();

	void run_1D(int num_iters = -1);
	void run_2D(int num_iters = -1);

	const RegularMesh * _mesh;
	SimulationData _simdata;

	// user defined data
	const double * _rel_permittivity;
	const double * _rel_permeability;
	const double * _current_density_x;
	const double * _field_multiplier_x;
	const double * _current_density_y;
	const double * _field_multiplier_y;
	const double * _current_density_z;
	const double * _field_multiplier_z;
	std::vector<SignalGenerator> _signals;
	std::vector<SignalGenerator> _modulators;


	// internal simulation variables
		bool _is_allocated, _include_current_density, _time_variable_coefficients;
		const double _eps0 = 8.854e-12;
		const double _mu0 = 4*3.14159e-7;
		const double _c0 = 2.99792458e+8;
		double _dx, _dt, _CourantFactor;
		double _nPML, _num_iters, _current_iter, _tcur;

		
		std::vector<double> _default_current_density_x;
		std::vector<double> _default_current_density_y;
		std::vector<double> _default_current_density_z;
		std::vector<double> _default_field_multiplier_x;
		std::vector<double> _default_field_multiplier_y;
		std::vector<double> _default_field_multiplier_z;
		

		// field stuff
		double * D_z;	// note that this is not the real D field, but it is normalized
		double * I_Hx;
		double * I_Hy;
		double * E_z;	// note that this is not the real E field, but it is normalized
		double * H_y;
		double * H_x;
		double * gaz;
		double * jnx; 
		double * jny;
		double * jnz;

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