#include "PlasmaDriftDiffusionSimulation.hpp"


using namespace std;

PlasmaDriftDiffusionSimulation::PlasmaDriftDiffusionSimulation(){
	
	_diffusion_coeff = nullptr;
	_gain_term = nullptr;
	_loss_term = nullptr;
	_init_density = nullptr;

	_mesh = nullptr;

	_dx = 0.0;
	_dt = 0.0;
	_CourantFactor = 0.5;
	_num_iters = 0;
	_current_iter = 0;
	_tcur = 0.0;

	_is_allocated = false;
}

PlasmaDriftDiffusionSimulation::~PlasmaDriftDiffusionSimulation(){

}

void PlasmaDriftDiffusionSimulation::bind_mesh(const RegularMesh & mesh){
	_mesh = &mesh;
	_dx = _mesh->res();
}

void PlasmaDriftDiffusionSimulation::bind_diffusion_coeff(const double * diffusion_coeff){
	_diffusion_coeff = diffusion_coeff;
}

void PlasmaDriftDiffusionSimulation::bind_gain_term(const double * gain_term){
	_gain_term = gain_term;
}

void PlasmaDriftDiffusionSimulation::bind_loss_term(const double * loss_term){
	_loss_term = loss_term;
}

void PlasmaDriftDiffusionSimulation::view_results(){
	// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_color_ramp(CRamp::MATLAB_PARULA);
		simvis.set_colorby_field("density");
		simvis.set_color_interpolation(false);
		simvis.set_frequency_Hz(30);
		simvis.run();
}

void PlasmaDriftDiffusionSimulation::output_HDF5(std::string outname){
	_simdata.write_HDF5(outname);
}

void PlasmaDriftDiffusionSimulation::run(int num_iters){
	
	if (!_is_allocated) {
		
		preRunCheck();
		allocate_fields();
		if (_dt == 0.0){
			//_dt = _CourantFactor * _dx * _dx / std::max_element(_diffusion_coeff, _diffusion_coeff + _mesh->nodecount());
			_dt = 0.5 * _CourantFactor * _dx * _dx / _diffusion_coeff[0];
		}
		allocate_simdata();
		_is_allocated = true;

		cout << "dx: " << _dx << endl;
		cout << "Courant Factor: " << _CourantFactor << endl;
		cout << "dt: " << _dt << endl;
	}
	//if (_is_allocated && _time_variable_coefficients) recalcCoeffs();

	
	if (_mesh->num_dims() == 1){
		cout << "1D Plasma Drift Diffusion is not yet plugged in!" << endl;
	}
	else if (_mesh->num_dims() == 2){
		run_2D(num_iters);
	}
	else if (_mesh->num_dims() == 3){
		cout << "3D Plasma Drift Diffusion is not yet implemented!" << endl;
	}
}

void PlasmaDriftDiffusionSimulation::run_2D(int num_iters){
	// discretized using a simple explicit euler scheme

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else {
		end_iter = _current_iter + num_iters;
		if (end_iter > _num_iters) end_iter = _num_iters;
	}
	
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;

		for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){
			for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){

				cind = _mesh->reg_inds_to_glob_ind(i, j);
				lind = _mesh->reg_inds_to_glob_ind(i-1, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);
				uind = _mesh->reg_inds_to_glob_ind(i, j+1);
				dind = _mesh->reg_inds_to_glob_ind(i, j-1);

				_density[cind] = _density_old[cind] + _dt*(_gain_term[cind]
												- _loss_term[cind]
												+ _diffusion_coeff[cind]/_dx/_dx*(_density_old[lind] + _density_old[rind]
													+ _density_old[uind] + _density_old[dind] - 4*_density_old[cind]));

			}
		}

		//cout << "density at center: " << _density.at(_mesh->nodecount()/2) << endl;

		// update the old density
		for (auto i=0; i<_mesh->nodecount(); i++) _density_old[i] = _density[i];
	
		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "density", _density.front());

		// update iteration count
		_current_iter++;

	}


}

void PlasmaDriftDiffusionSimulation::preRunCheck(){

}


void PlasmaDriftDiffusionSimulation::allocate_fields(){

	_density.assign(_mesh->nodecount(), 0);
	_density_old.assign(_mesh->nodecount(), 0);

	if (_init_density != nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++){
			_density.at(i) = _init_density[i];
			_density_old.at(i) = _init_density[i];
		}
	}

	// diffusion coefficient
	if (_diffusion_coeff == nullptr){
		_default_diffusion_coeff.assign(_mesh->nodecount(), (_k_b*2)/(_m_e*4.4e+12)); // (_k_b*3.20435313e-19)/(_m_e*4.4e+12)
		_diffusion_coeff = &_default_diffusion_coeff.front();
	}

	// gain term
	if (_gain_term == nullptr){
		_default_gain_term.assign(_mesh->nodecount(), 0.0);
		_gain_term = &_default_gain_term.front();
	}

	// loss term
	if (_loss_term == nullptr){
		_default_loss_term.assign(_mesh->nodecount(), 0.0);
		_loss_term = &_default_loss_term.front();
	}
}

void PlasmaDriftDiffusionSimulation::allocate_simdata(){
	// set up the simulation data holder
	_simdata.bind_mesh(*_mesh);
	_simdata.add_field("density");
	_simdata.set_time_span(0.0, _dt, _num_iters*_dt);
	_simdata.print_summary();
}
