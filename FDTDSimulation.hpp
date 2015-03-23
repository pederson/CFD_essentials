#ifndef _FDTDSIMULATION_H
#define _FDTDSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerSimulation.hpp"
#include "SignalGenerator.hpp"
#include "cvector.hpp"

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
	const double & E_z() const {return _E_z[0];};
	const double & H_y() const {return _H_y[0];};
	const double & H_x() const {return _H_x[0];};
	const double & E_z(double t) {return _simdata.get_data_at_time(t, "E_z");}; 
	const double * E_z_ptr() const {return _E_z;};
	const double * H_x_ptr() const {return _H_x;};
	const double * H_y_ptr() const {return _H_y;};
	

	// mutators
	void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void add_gaussian_source(double t0, double spread, double xloc, double yloc=0, double zloc=0);
	void add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc=0, double zloc=0);
	void add_gaussian_modulator(unsigned int signal_idx, double t0, double spread);
	void add_sinusoidal_modulator(unsigned int signal_idx, double freq_Hz, double phase=0.0);
	
	void bind_mesh(const RegularMesh & mesh);
	//void bind_metal_nodes(const double * metal_nodes);
	void bind_rel_permittivity(const double * rel_permittivity);
	void bind_rel_permittivity(const cvector & rel_permittivity_cv);
	//void bind_rel_permeability(const double * rel_permeability);
	void bind_current_density_x(const double * current_density_x);
	void bind_current_density_y(const double * current_density_y);
	void bind_current_density_z(const double * current_density_z);

	void set_num_iters(unsigned int num_iters) {_num_iters = num_iters;};
	//void set_invariant_time(unsigned int timespan);
	//void set_invariant_resolution(unsigned int res);


	void view_results();
	void output_HDF5(std::string outname="FDTD_output.h5");

	void run(int num_iters = -1);

private:

	void preRunCheck();
	void prepareModulator();
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
	const double * _current_density_y;
	const double * _current_density_z;
	cvector _rel_permittivity_cv;
	cvector _rel_permeability_cv;
	cvector _current_density_x_cv;
	cvector _current_density_y_cv;
	cvector _current_density_z_cv;
	std::vector<SignalGenerator> _signals;
	std::vector<SignalGenerator> _modulators;


	// internal simulation variables
		bool _is_allocated, _include_current_density, _time_variable_coefficients;
		const double _pi = 3.14159265358979323846264338327950288;
		const double _eps0 = 8.854e-12;
		const double _mu0 = _pi*4.0e-7;
		const double _c0 = 2.99792458e+8;
		double _dx, _dt, _CourantFactor;
		double _nPML, _num_iters, _current_iter, _tcur;

		
		std::vector<double> _default_rel_permittivity;
		std::vector<double> _default_rel_permeability;
		std::vector<double> _default_current_density_x;
		std::vector<double> _default_current_density_y;
		std::vector<double> _default_current_density_z;


		SignalGenerator _source_modulator; // modulates all sources with tanh
		double _source_modulator_width; 	// how many time steps until it reaches max value
		

		// field stuff
		double * _Dn_z;	// note that this is not the real D field, but it is normalized
		double * _I_Hx;
		double * _I_Hy;
		double * _En_z;	// note that this is not the real E field, but it is normalized
		double * _E_z;	// this is the real E field
		double * _H_y;
		double * _H_x;
		//double * _gaz;
		//double * _jnx; 
		//double * _jny;
		//double * _jnz;

		// PML stuff only used internally
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