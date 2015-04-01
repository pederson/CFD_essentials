#ifndef _FDTDSIMULATION_H
#define _FDTDSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerSimulation.hpp"
#include "SignalGenerator.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <functional>

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
	const SimulationData & simdata() const {return _simdata;};
	const double dt() const {return _dt;};

	// mutators
	void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void add_gaussian_source(double t0, double spread, double xloc, double yloc=0, double zloc=0);
	void add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc=0, double zloc=0);
	void add_gaussian_modulator(unsigned int signal_idx, double t0, double spread);
	void add_sinusoidal_modulator(unsigned int signal_idx, double freq_Hz, double phase=0.0);
	void bind_mesh(const RegularMesh & mesh);
	//void bind_metal_nodes(const double * metal_nodes);
	void bind_rel_permittivity(const double * rel_permittivity);
	void bind_conductivity(const double * conductivity);
	void bind_single_pole(const double * numerator, const double * frequency_pole);
	//void bind_rel_permeability(const double * rel_permeability);
	void bind_current_density_x(const double * current_density_x);
	void bind_current_density_y(const double * current_density_y);
	void bind_current_density_z(const double * current_density_z);
	//void bind_metal_nodes(std::function<double(unsigned int) metal_nodes_fn);
	void bind_rel_permittivity(std::function<double(unsigned int)> rel_permittivity_fn);
	void bind_conductivity(std::function<double(unsigned int)> conductivity_fn);
	void bind_single_pole(std::function<double(unsigned int)> numerator_fn, std::function<double(unsigned int)> frequency_pole_fn);
	//void bind_rel_permeability(std::function<double(unsigned int) rel_permeability_fn);
	void bind_current_density_x(std::function<double(unsigned int)> current_density_x_fn);
	void bind_current_density_y(std::function<double(unsigned int)> current_density_y_fn);
	void bind_current_density_z(std::function<double(unsigned int)> current_density_z_fn);

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
	void allocate_coeffs();
	void allocate_PML();
	void allocate_simdata();

	void run_1D(int num_iters = -1);
	void run_2D(int num_iters = -1);
	void run_3D(int num_iters = -1);

	const RegularMesh * _mesh;
	SimulationData _simdata;

	// user defined data
	const double * _rel_permittivity;
	const double * _rel_permeability;
	const double * _conductivity;
	const double * _permittivity_single_pole_freq;
	const double * _permittivity_single_pole_numerator;
	const double * _current_density_x;
	const double * _current_density_y;
	const double * _current_density_z;

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

		// functions used for defining properties
		std::function<double(unsigned int)> _rel_permittivity_fn;
		std::function<double(unsigned int)> _conductivity_fn;
		std::function<double(unsigned int)> _permittivity_single_pole_freq_fn;
		std::function<double(unsigned int)> _permittivity_single_pole_numerator_fn;
		std::function<double(unsigned int)> _rel_permeability_fn;
		std::function<double(unsigned int)> _current_density_x_fn;
		std::function<double(unsigned int)> _current_density_y_fn;
		std::function<double(unsigned int)> _current_density_z_fn;

		SignalGenerator _source_modulator; // modulates all sources with tanh
		double _source_modulator_width; 	// how many time steps until it reaches max value
		

		// field stuff
		// Electric polarization field
		double * _Dn_x;	// note that this is the normalized D field
		double * _Dn_y;
		double * _Dn_z;	
		
		// integrated quantities
		double * _I_Hx;
		double * _I_Hy;
		double * _I_Hz;
		double * _I_Dnx;
		double * _I_Dny;
		double * _I_Dnz;

		// Electric fields
		double * _En_x;	// note that this is the normalized E field
		double * _En_y;
		double * _En_z;	
		double * _E_x;	// real E fields
		double * _E_y;
		double * _E_z;
		
		// Magnetic fields
		double * _H_x;	// real H fields
		double * _H_y;
		double * _H_z;

		// used for frequency dependent materials
		double * _I_Ex;	// non-zero for non-zero conductivity
		double * _I_Ey;
		double * _I_Ez;
		double * _S_n_x;	// non-zero for non-zero single pole permittivity
		double * _S_n_y;
		double * _S_n_z;
		double * _S_nm1_x;// non-zero for non-zero single pole permittivity
		double * _S_nm1_y;
		double * _S_nm1_z;

		// PML stuff only used internally
		double * gi1;
		double * gi2;
		double * gi3;
		double * gj1;
		double * gj2;
		double * gj3;
		double * gk1;
		double * gk2;
		double * gk3;
		double * fi1;
		double * fi2;
		double * fi3;
		double * fj1;
		double * fj2;
		double * fj3;
		double * fk1;
		double * fk2;
		double * fk3;



};





#endif

