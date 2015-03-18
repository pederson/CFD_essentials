#ifndef _PLASMAMOMENTUMSIMULATION_H
#define _PLASMAMOMENTUMSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerSimulation.hpp"
#include "SignalGenerator.hpp"

#include <vector>
#include <string>
#include <iostream>

class PlasmaMomentumSimulation{
public:

	PlasmaMomentumSimulation();
	~PlasmaMomentumSimulation();

	// inspectors
	const double & velocity_x() const {return _velocity_x.front();};
	const double & velocity_y() const {return _velocity_y.front();};
	const double & velocity_z() const {return _velocity_z.front();};
	const double & collision_rate() const {return _collision_rate[0];};
	const double & E_x() const {return _E_x[0];};
	const double & E_y() const {return _E_y[0];};
	const double & E_z() const {return _E_z[0];};

	// mutators
	//void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
	void set_initial_velocity_x(const double * init_vel_x) {_init_vel_x = init_vel_x;};
	void set_initial_velocity_y(const double * init_vel_y) {_init_vel_y = init_vel_y;};
	void set_initial_velocity_z(const double * init_vel_z) {_init_vel_z = init_vel_z;};
	void set_time_step(double dt) {_dt = dt;};
	void bind_mesh(const RegularMesh & mesh);
	void bind_collision_rate(const double * diffusion_coeff);
	void bind_E_x(const double * E_x);
	void bind_E_y(const double * E_y);
	void bind_E_z(const double * E_z);
	void bind_E_x_prev(const double * E_x_prev);
	void bind_E_y_prev(const double * E_y_prev);
	void bind_E_z_prev(const double * E_z_prev);
	void set_num_iters(unsigned int num_iters) {_num_iters = num_iters;};
	
	// other
	void view_results();
	void output_HDF5(std::string outname="PlasmaMomentum_output.h5");

	void run(int num_iters = -1);

private:

	void preRunCheck();
	void allocate_fields();
	void allocate_simdata();

	void run_2D(int num_iters = -1);

	const RegularMesh * _mesh;
	SimulationData _simdata;

	// user defined data
	const double * _collision_rate;	// if it varies over the field
	const double * _E_x, * _E_x_prev;
	const double * _E_y, * _E_y_prev;
	const double * _E_z, * _E_z_prev;

	const double * _init_vel_x;
	const double * _init_vel_y;
	const double * _init_vel_z;


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

		// density field
		std::vector<double> _velocity_x;
		std::vector<double> _velocity_y;
		std::vector<double> _velocity_z;

		// defaults
		std::vector<double> _default_collision_rate;
		std::vector<double> _default_E_x;
		std::vector<double> _default_E_y;
		std::vector<double> _default_E_z;

};


#endif