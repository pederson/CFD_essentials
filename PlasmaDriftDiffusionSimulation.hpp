#ifndef _PLASMADRIFTDIFFUSIONSIMULATION_H
#define _PLASMADRIFTDIFFUSIONSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerSimulation.hpp"
#include "SignalGenerator.hpp"

#include <vector>
#include <string>
#include <iostream>

class PlasmaDriftDiffusionSimulation{
public:

	PlasmaDriftDiffusionSimulation();
	~PlasmaDriftDiffusionSimulation();

	// inspectors
	const double & density() const {return _density.front();};
	const double * density_ptr() const {return &_density.front();};
	//const double & diffusion_coeff() const {return _diffusion_coeff[0];};
	//const double & gain_term() const {return _gain_term[0];};
	//const double & loss_term() const {return _loss_term[0];};

	// mutators
	//void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void set_initial_density(const double * init_density);
	void set_time_step(double dt) {_dt = dt;};
	void bind_mesh(const RegularMesh & mesh);
	void bind_diffusion_coeff(const double * diffusion_coeff);
	void bind_gain_term(const double * gain_term);
	void bind_loss_term(const double * loss_term);
	void set_num_iters(unsigned int num_iters) {_num_iters = num_iters;};
	
	// other
	void view_results();
	void output_HDF5(std::string outname="PlasmaDriftDiffusion_output.h5");

	void run(int num_iters = -1);

private:

	void preRunCheck();
	void allocate_fields();
	void allocate_simdata();

	void run_2D(int num_iters = -1);

	const RegularMesh * _mesh;
	SimulationData _simdata;

	// user defined data
	const double * _diffusion_coeff;	// if it varies over the field
	const double * _gain_term;
	const double * _loss_term;
	const double * _init_density;
	//const cvector _diffusion_coeff_cv, _gain_term_cv, _loss_term_cv;


	// internal simulation variables
		bool _is_allocated;
		const double _pi = 3.14159265358979323846264338327950288;
		const double _eps0 = 8.854e-12;
		const double _mu0 = _pi*4.0e-7;

		const double _m_e = 9.10938291e-31;		// electron mass
		const double _q_e = 1.60217657e-19;		// electron charge
		const double _k_b = 1.3806488e-23;		// Boltzmann constant

		const double _c0 = 2.99792458e+8;
		double _dx, _dt, _CourantFactor;
		double _num_iters, _current_iter, _tcur;

		
		SignalGenerator _source_modulator; // modulates all sources with tanh
		double _source_modulator_width; 	// how many time steps until it reaches max value
		
		// density field
		std::vector<double> _density, _density_old;	

		// defaults
		std::vector<double> _default_diffusion_coeff;
		std::vector<double> _default_gain_term;
		std::vector<double> _default_loss_term;

};


#endif
