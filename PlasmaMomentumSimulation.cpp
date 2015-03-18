#include "PlasmaMomentumSimulation.hpp"


using namespace std;

PlasmaMomentumSimulation::PlasmaMomentumSimulation(){
	
	_collision_rate = nullptr;
	_E_x = nullptr;
	_E_y = nullptr;
	_E_z = nullptr;
	_E_x_prev = nullptr;
	_E_y_prev = nullptr;
	_E_z_prev = nullptr;

	_init_vel_x = nullptr;
	_init_vel_y = nullptr;
	_init_vel_z = nullptr;

	_mesh = nullptr;

	_dx = 0.0;
	_dt = 0.0;
	_CourantFactor = 0.5;
	_num_iters = 0;
	_current_iter = 0;
	_tcur = 0.0;

	_is_allocated = false;
}

PlasmaMomentumSimulation::~PlasmaMomentumSimulation(){

}

void PlasmaMomentumSimulation::bind_mesh(const RegularMesh & mesh){
	_mesh = &mesh;
	_dx = _mesh->res();
}

void PlasmaMomentumSimulation::bind_collision_rate(const double * collision_rate){
	_collision_rate = collision_rate;
}

void PlasmaMomentumSimulation::bind_E_x(const double * E_x){
	_E_x = E_x;
}

void PlasmaMomentumSimulation::bind_E_y(const double * E_y){
	_E_y = E_y;
}

void PlasmaMomentumSimulation::bind_E_z(const double * E_z){
	_E_z = E_z;
}

void PlasmaMomentumSimulation::bind_E_x_prev(const double * E_x_prev){
	_E_x_prev = E_x_prev;
}

void PlasmaMomentumSimulation::bind_E_y_prev(const double * E_y_prev){
	_E_y_prev = E_y_prev;
}

void PlasmaMomentumSimulation::bind_E_z_prev(const double * E_z_prev){
	_E_z_prev = E_z_prev;
}

void PlasmaMomentumSimulation::view_results(){
	// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_color_ramp(CRamp::MATLAB_PARULA);
		simvis.set_colorby_field("velocity_x");
		simvis.set_color_interpolation(false);
		simvis.set_frequency_Hz(30);
		simvis.run();
}

void PlasmaMomentumSimulation::output_HDF5(std::string outname){
	_simdata.write_HDF5(outname);
}

void PlasmaMomentumSimulation::run(int num_iters){
	
	if (!_is_allocated) {
		
		preRunCheck();
		allocate_fields();
		if (_dt == 0.0){
			_dt = 0.5 * _CourantFactor * _dx;
		}
		allocate_simdata();
		_is_allocated = true;

		cout << "dx: " << _dx << endl;
		cout << "Courant Factor: " << _CourantFactor << endl;
		cout << "dt: " << _dt << endl;
	}
	//if (_is_allocated && _time_variable_coefficients) recalcCoeffs();

	
	if (_mesh->num_dims() == 1){
		cout << "1D Plasma Momentum is not yet plugged in!" << endl;
	}
	else if (_mesh->num_dims() == 2){
		run_2D(num_iters);
	}
	else if (_mesh->num_dims() == 3){
		cout << "3D Plasma Momentum is not yet implemented!" << endl;
	}
}

void PlasmaMomentumSimulation::run_2D(int num_iters){
	// discretized using Crank Nicholson scheme

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else {
		end_iter = _current_iter + num_iters;
		if (end_iter > _num_iters) end_iter = _num_iters;
	}
	
	unsigned int cind;
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;

		for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){
			for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){

				cind = _mesh->reg_inds_to_glob_ind(i, j);

				_velocity_x[cind] = _velocity_x[cind]*(1-_CourantFactor*_dt*_collision_rate[cind])/(1+_CourantFactor*_dt*_collision_rate[cind]) 
									- (_E_x[cind] + _E_x_prev[cind])*0.5*(_q_e*_dt)/(_m_e*(1+_CourantFactor*_dt*_collision_rate[cind]));

				_velocity_y[cind] = _velocity_y[cind]*(1-_CourantFactor*_dt*_collision_rate[cind])/(1+_CourantFactor*_dt*_collision_rate[cind]) 
									- (_E_y[cind] + _E_y_prev[cind])*0.5*(_q_e*_dt)/(_m_e*(1+_CourantFactor*_dt*_collision_rate[cind]));

			}
		}

		//cout << "density at center: " << _density.at(_mesh->nodecount()/2) << endl;

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "velocity_x", _velocity_x.front());
		_simdata.add_data_at_index(n, "velocity_y", _velocity_y.front());

		// update iteration count
		_current_iter++;

	}


}

void PlasmaMomentumSimulation::preRunCheck(){

}


void PlasmaMomentumSimulation::allocate_fields(){

	// x velocity
	_velocity_x.assign(_mesh->nodecount(), 0);
	if (_init_vel_x != nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++){
			_velocity_x.at(i) = _init_vel_x[i];
		}
	}

	// collision rate
	if (_collision_rate == nullptr){
		_default_collision_rate.assign(_mesh->nodecount(), 0); 
		_collision_rate = &_default_collision_rate.front();
	}

	// x Efield
	if (_E_x == nullptr){
		_default_E_x.assign(_mesh->nodecount(), 0.0);
		_E_x = &_default_E_x.front();
	}

	// x Efield previous
	if (_E_x_prev == nullptr){
		_E_x_prev = _E_x;
	}

	if (_mesh->num_dims() > 1){
		// y velocity
		_velocity_y.assign(_mesh->nodecount(), 0);
		if (_init_vel_y != nullptr){
			for (auto i=0; i<_mesh->nodecount(); i++){
				_velocity_y.at(i) = _init_vel_y[i];
			}
		}

		// y Efield
		if (_E_y == nullptr){
			_default_E_y.assign(_mesh->nodecount(), 0.0);
			_E_y = &_default_E_y.front();
		}

		// y Efield previous
		if (_E_y_prev == nullptr){
			_E_y_prev = _E_y;
		}
	}

	if (_mesh->num_dims() > 2){
		// z velocity
		_velocity_z.assign(_mesh->nodecount(), 0);
		if (_init_vel_z != nullptr){
			for (auto i=0; i<_mesh->nodecount(); i++){
				_velocity_z.at(i) = _init_vel_z[i];
			}
		}

		// z Efield
		if (_E_z == nullptr){
			_default_E_z.assign(_mesh->nodecount(), 0.0);
			_E_z = &_default_E_z.front();
		}

		// z Efield previous
		if (_E_z_prev == nullptr){
			_E_z_prev = _E_z;
		}
	}

}

void PlasmaMomentumSimulation::allocate_simdata(){
	// set up the simulation data holder
	_simdata.bind_mesh(*_mesh);
	_simdata.add_field("velocity_x");
	_simdata.add_field("velocity_y");
	_simdata.set_time_span(0.0, _dt, _num_iters*_dt);
	_simdata.print_summary();
}