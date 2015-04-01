#include "FDTDSimulation.hpp"

using namespace std;

FDTDSimulation::FDTDSimulation(){
	_Dn_x = nullptr;
	_Dn_y = nullptr;
	_Dn_z = nullptr;

	_I_Hx = nullptr;
	_I_Hy = nullptr;
	_I_Hz = nullptr;
	_I_Dnx = nullptr;
	_I_Dny = nullptr;
	_I_Dnz = nullptr;

	_En_x = nullptr;
	_En_y = nullptr;
	_En_z = nullptr;
	_E_x = nullptr;
	_E_y = nullptr;
	_E_z = nullptr;

	_H_x = nullptr;
	_H_y = nullptr;
	_H_z = nullptr;

	_I_Ex = nullptr;
	_I_Ey = nullptr;
	_I_Ez = nullptr;
	_S_n_x = nullptr;
	_S_n_y = nullptr;
	_S_n_z = nullptr;
	_S_nm1_x = nullptr;
	_S_nm1_y = nullptr;
	_S_nm1_z = nullptr;

	gi1 = nullptr;
	gi2 = nullptr;
	gi3 = nullptr;
	gj1 = nullptr;
	gj2 = nullptr;
	gj3 = nullptr;
	gk1 = nullptr;
	gk2 = nullptr;
	gk3 = nullptr;
	fi1 = nullptr;
	fi2 = nullptr;
	fi3 = nullptr;
	fj1 = nullptr;
	fj2 = nullptr;
	fj3 = nullptr;
	fk1 = nullptr;
	fk2 = nullptr;
	fk3 = nullptr;

	_rel_permittivity = nullptr;
	_rel_permeability = nullptr;
	_conductivity = nullptr;
	_permittivity_single_pole_freq = nullptr;
	_permittivity_single_pole_numerator = nullptr;
	_current_density_x = nullptr;
	_current_density_y = nullptr;
	_current_density_z = nullptr;

	_rel_permittivity_fn = nullptr;
	_conductivity_fn = nullptr;
	_permittivity_single_pole_freq_fn = nullptr;
	_permittivity_single_pole_numerator_fn = nullptr;
	_rel_permeability_fn = nullptr;
	_current_density_x_fn = nullptr;
	_current_density_y_fn = nullptr;
	_current_density_z_fn = nullptr;

	_mesh = nullptr;

	_dx = 0.0;
	_dt = 0.0;
	_CourantFactor = 0.5;
	_nPML = 8;
	_num_iters = 0;
	_current_iter = 0;
	_tcur = 0.0;

	_is_allocated = false;
	_include_current_density = false;
	_time_variable_coefficients = false;


}

FDTDSimulation::~FDTDSimulation(){
	if (_Dn_x != nullptr) delete[] _Dn_x;
	if (_Dn_y != nullptr) delete[] _Dn_y;
	if (_Dn_z != nullptr) delete[] _Dn_z;

	if (_I_Hx != nullptr) delete[] _I_Hx;
	if (_I_Hy != nullptr) delete[] _I_Hy;
	if (_I_Hz != nullptr) delete[] _I_Hz;
	if (_I_Dnx != nullptr) delete[] _I_Dnx;
	if (_I_Dny != nullptr) delete[] _I_Dny;
	if (_I_Dnz != nullptr) delete[] _I_Dnz;

	if (_En_x != nullptr) delete[] _En_x;
	if (_En_y != nullptr) delete[] _En_y;
	if (_En_z != nullptr) delete[] _En_z;
	if (_E_x != nullptr) delete[] _E_x;
	if (_E_y != nullptr) delete[] _E_y;
	if (_E_z != nullptr) delete[] _E_z;
	
	if (_H_x != nullptr) delete[] _H_x;
	if (_H_y != nullptr) delete[] _H_y;
	if (_H_z != nullptr) delete[] _H_z;

	if (_I_Ex != nullptr) delete[] _I_Ex;
	if (_I_Ey != nullptr) delete[] _I_Ey;
	if (_I_Ez != nullptr) delete[] _I_Ez;
	if (_S_n_x != nullptr) delete[] _S_n_x;
	if (_S_n_y != nullptr) delete[] _S_n_y;
	if (_S_n_z != nullptr) delete[] _S_n_z;
	if (_S_nm1_x != nullptr) delete[] _S_nm1_x;
	if (_S_nm1_y != nullptr) delete[] _S_nm1_y;
	if (_S_nm1_z != nullptr) delete[] _S_nm1_z;

	if (gi1 != nullptr) delete[] gi1;
	if (gi2 != nullptr) delete[] gi2;
	if (gi3 != nullptr) delete[] gi3;
	if (gj1 != nullptr) delete[] gj1;
	if (gj2 != nullptr) delete[] gj2;
	if (gj3 != nullptr) delete[] gj3;
	if (gk1 != nullptr) delete[] gk1;
	if (gk2 != nullptr) delete[] gk2;
	if (gk3 != nullptr) delete[] gk3;
	if (fi1 != nullptr) delete[] fi1;
	if (fi2 != nullptr) delete[] fi2;
	if (fi3 != nullptr) delete[] fi3;
	if (fj1 != nullptr) delete[] fj1;
	if (fj2 != nullptr) delete[] fj2;
	if (fj3 != nullptr) delete[] fj3;
	if (fk1 != nullptr) delete[] fk1;
	if (fk2 != nullptr) delete[] fk2;
	if (fk3 != nullptr) delete[] fk3;
}
void FDTDSimulation::set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers){

}

void FDTDSimulation::add_gaussian_source(double t0, double spread, double xloc, double yloc, double zloc){
	SignalGenerator sig, mod;
	sig.set_gaussian(spread, t0);
	sig.set_location(xloc, yloc, zloc);
	mod.set_location(xloc, yloc, zloc);

	_signals.push_back(sig);
	_modulators.push_back(mod);

}

void FDTDSimulation::add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc, double zloc){
	SignalGenerator sig, mod;
	sig.set_sinusoid(freq_Hz, phase);
	sig.set_location(xloc, yloc, zloc);
	mod.set_location(xloc, yloc, zloc);

	_signals.push_back(sig);
	_modulators.push_back(mod);
}

void FDTDSimulation::add_gaussian_modulator(unsigned int signal_idx, double t0, double spread){
	SignalGenerator mod;
	mod.set_gaussian(spread, t0);

	_modulators.at(signal_idx) = mod;
}

void FDTDSimulation::add_sinusoidal_modulator(unsigned int signal_idx, double freq_Hz, double phase){
	SignalGenerator mod;
	mod.set_sinusoid(freq_Hz, phase);

	_modulators.at(signal_idx) = mod;
}

void FDTDSimulation::bind_mesh(const RegularMesh & mesh) {
	_mesh = &mesh;
	_dx = _mesh->res();
	_dt = _CourantFactor*_dx/_c0;
}

void FDTDSimulation::bind_rel_permittivity(const double * rel_permittivity){
	_rel_permittivity = rel_permittivity;
}

void FDTDSimulation::bind_conductivity(const double * conductivity){
	_conductivity = conductivity;
}

void FDTDSimulation::bind_single_pole(const double * numerator, const double * freq_pole){
	_permittivity_single_pole_freq = freq_pole;
	_permittivity_single_pole_numerator = numerator;
}

void FDTDSimulation::bind_current_density_x(const double * current_density_x){
	_current_density_x = current_density_x;
}

void FDTDSimulation::bind_current_density_y(const double * current_density_y){
	_current_density_y = current_density_y;
}

void FDTDSimulation::bind_current_density_z(const double * current_density_z){
	_current_density_z = current_density_z;
}

void FDTDSimulation::bind_rel_permittivity(function<double(unsigned int)> rel_permittivity_fn){
	_rel_permittivity_fn = rel_permittivity_fn;
}

void FDTDSimulation::bind_conductivity(function<double(unsigned int)> conductivity_fn){
	_conductivity_fn = conductivity_fn;
}

void FDTDSimulation::bind_single_pole(function<double(unsigned int)> numerator_fn, function<double(unsigned int)> frequency_pole_fn){
	_permittivity_single_pole_numerator_fn = numerator_fn;
	_permittivity_single_pole_freq_fn = frequency_pole_fn;
}

void FDTDSimulation::bind_current_density_x(function<double(unsigned int)> current_density_x_fn){
	_current_density_x_fn = current_density_x_fn;
}

void FDTDSimulation::bind_current_density_y(function<double(unsigned int)> current_density_y_fn){
	_current_density_y_fn = current_density_y_fn;
}

void FDTDSimulation::bind_current_density_z(function<double(unsigned int)> current_density_z_fn){
	_current_density_z_fn = current_density_z_fn;
}


void FDTDSimulation::view_results(){

		// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_color_ramp(CRamp::MATLAB_PARULA);
		simvis.set_colorby_field("E_z");
		//simvis.set_colorby_field("H_y");
		simvis.set_color_alpha(_rel_permittivity);
		simvis.set_color_interpolation(false);
		simvis.set_snapshot_increment(3);
		simvis.run();

}

void FDTDSimulation::output_HDF5(std::string outname){
	_simdata.write_HDF5(outname);
}

void FDTDSimulation::run(int num_iters){

	if (!_is_allocated) {
		preRunCheck();
		allocate_fields();
		allocate_coeffs();
		allocate_PML();
		allocate_simdata();
		_is_allocated = true;

		_dt = _CourantFactor*_dx/_c0;

		cout << "dx: " << _dx << endl;
		cout << "Courant Factor: " << _CourantFactor << endl;
		cout << "dt: " << _dt << endl;

		cout << "starting the run in " << _mesh->num_dims() << "D" << endl;
	}

	if (_mesh->num_dims() == 1){
		cout << "1D FDTD does not support complex materials yet!" << endl;
		run_1D(num_iters);
	}
	else if (_mesh->num_dims() == 2){
		run_2D(num_iters);
	}
	else if (_mesh->num_dims() == 3){
		cout << "3D FDTD is not yet tested!" << endl;
		run_3D(num_iters);
	}

}

void FDTDSimulation::run_1D(int num_iters){

	
	double pulse;
	double _srcx;

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else {
		end_iter = _current_iter + num_iters;
		if (end_iter > _num_iters) end_iter = _num_iters;
	}
	
	double curle, sourcemodval;
	unsigned int cind, lind, rind;
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;
		sourcemodval = _source_modulator.value(double(_current_iter)/_source_modulator_width*SIGNALGENERATOR_PI/2.0);
		//cout << endl;

		// update D field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols

			cind = _mesh->reg_inds_to_glob_ind(i);
			lind = _mesh->reg_inds_to_glob_ind(i-1);
			rind = _mesh->reg_inds_to_glob_ind(i+1);

			
			_Dn_z[cind] = gi3[i]*_Dn_z[cind]
					  + gi2[i]*_CourantFactor*(_H_y[cind] - _H_y[lind])
					  - _CourantFactor*_dx*_current_density_z_fn(cind)*sourcemodval;
			//_Dn_z[cind] = _Dn_z[cind]
			//		  + _CourantFactor*(_H_y[cind] - _H_y[lind]);

		}

		//cout << "about to impose source" << endl;
		// impose sources
		for (auto i=0; i<_signals.size(); i++){
			pulse = _signals.at(i).value(_tcur)*_modulators.at(i).value(_tcur) * sourcemodval;
			_srcx = _signals.at(i).xloc();
			_Dn_z[_mesh->nearest_node(_srcx)] = pulse;
			//cout << "source value: " << pulse << endl;
		}
		//cout << "imposed source" << endl;
		
		// calculate Ez field
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){ // cols
			cind = _mesh->reg_inds_to_glob_ind(i);

			// with pml, conductivity, and single pole
			_En_z[cind] = (_Dn_z[cind] - _I_Ez[cind] - exp(-_dt*_permittivity_single_pole_freq_fn(cind))*_S_nm1_z[cind])/
						  (_rel_permittivity_fn(cind) + _conductivity_fn(cind)*_dt/_eps0 + _permittivity_single_pole_numerator_fn(cind)*_dt);

			
			_I_Ez[cind] = _I_Ez[cind] + _conductivity_fn(cind)*_dt/_eps0*_En_z[cind];
			_S_n_z[cind] = exp(-_permittivity_single_pole_freq_fn(cind)*_dt)*_S_nm1_z[cind] + _permittivity_single_pole_numerator_fn(cind)*_dt*_En_z[cind];
			_S_nm1_z[cind] = _S_n_z[cind];

			_E_z[cind] = _En_z[cind] * sqrt(_mu0/_eps0); // this can be made faster
		}

		// need generalized boundary conditions here
		// set Ez at edges to zero for pml
		lind = _mesh->reg_inds_to_glob_ind(0);
		rind = _mesh->reg_inds_to_glob_ind(_mesh->reg_num_nodes_x()-1);
		_En_z[lind] = 0.0;
		_En_z[rind] = 0.0;
		

		// Calculate Hy
		for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
			cind = _mesh->reg_inds_to_glob_ind(i);
			rind = _mesh->reg_inds_to_glob_ind(i+1);

			curle = _En_z[rind] - _En_z[cind];
			_I_Hy[cind] = _I_Hy[cind] + curle;
			_H_y[cind] = fi3[i]*_H_y[cind]
					  + fi2[i]*_CourantFactor*curle;// + fi1[i]*_I_Hy[cind];
			//_H_y[cind] = _H_y[cind]
			//		  + _CourantFactor*(_En_z[rind] - _En_z[cind]);

		}
		

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", _E_z[0]);
		//_simdata.add_data_at_index(n, "H_y", _H_y[0]);

		// update iteration count
		_current_iter++;

	}


}

void FDTDSimulation::run_2D(int num_iters){

	
	double pulse;
	double _srcx, _srcy;

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else {
		end_iter = _current_iter + num_iters;
		if (end_iter > _num_iters) end_iter = _num_iters;
	}
	
	double curle, sourcemodval;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;
		sourcemodval = _source_modulator.value(double(_current_iter)/_source_modulator_width*SIGNALGENERATOR_PI/2.0);
		//cout << endl;

		// update D field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){

				cind = _mesh->reg_inds_to_glob_ind(i, j);
				lind = _mesh->reg_inds_to_glob_ind(i-1, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);
				uind = _mesh->reg_inds_to_glob_ind(i, j+1);
				dind = _mesh->reg_inds_to_glob_ind(i, j-1);

				
				_Dn_z[cind] = gi3[i]*gj3[j]*_Dn_z[cind]
						  + gi2[i]*gj2[j]*_CourantFactor*(_H_y[cind] - _H_y[lind] - _H_x[cind] + _H_x[dind])
						  - _CourantFactor*_dx*_current_density_z_fn(cind)*sourcemodval;

				
			}
		}

		//cout << "about to impose source" << endl;
		// impose sources
		for (auto i=0; i<_signals.size(); i++){
			pulse = _signals.at(i).value(_tcur)*_modulators.at(i).value(_tcur) * sourcemodval;
			_srcx = _signals.at(i).xloc();
			_srcy = _signals.at(i).yloc();
			_Dn_z[_mesh->nearest_node(_srcx, _srcy)] = pulse;
		}
		//cout << "imposed source" << endl;
		
		// calculate Ez field
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);

				// with pml, conductivity, and single pole
				_En_z[cind] = (_Dn_z[cind] - _I_Ez[cind] - exp(-_dt*_permittivity_single_pole_freq_fn(cind))*_S_nm1_z[cind])/
							  (_rel_permittivity_fn(cind) + _conductivity_fn(cind)*_dt/_eps0 + _permittivity_single_pole_numerator_fn(cind)*_dt);

				
				_I_Ez[cind] = _I_Ez[cind] + _conductivity_fn(cind)*_dt/_eps0*_En_z[cind];
				_S_n_z[cind] = exp(-_permittivity_single_pole_freq_fn(cind)*_dt)*_S_nm1_z[cind] + _permittivity_single_pole_numerator_fn(cind)*_dt*_En_z[cind];
				_S_nm1_z[cind] = _S_n_z[cind];

				_E_z[cind] = _En_z[cind] * sqrt(_mu0/_eps0); // this can be made faster
			
			}
		}

		// need generalized boundary conditions here
		// set Ez at edges to zero for pml
		for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
			lind = _mesh->reg_inds_to_glob_ind(0, j);
			rind = _mesh->reg_inds_to_glob_ind(_mesh->reg_num_nodes_x()-1, j);
			_En_z[lind] = 0.0;
			_En_z[rind] = 0.0;
		}
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			uind = _mesh->reg_inds_to_glob_ind(i, _mesh->reg_num_nodes_y()-1);
			dind = _mesh->reg_inds_to_glob_ind(i, 0);
			_En_z[uind] = 0.0;
			_En_z[dind] = 0.0;
		}
		

		// Calculate Hx
		for (auto j=0; j<_mesh->reg_num_nodes_y()-1; j++){
			for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);
				uind = _mesh->reg_inds_to_glob_ind(i, j+1);

				curle = _En_z[cind] - _En_z[uind];
				/*
				_I_Hx[cind] = _I_Hx[cind] + fi1[i]*curle;
				_H_x[cind] = fj3[j]*_H_x[cind]
						  + fj2[j]*_CourantFactor*(curle + _I_Hx[cind]);
				*/
				
				_I_Hx[cind] = _I_Hx[cind] + curle;
				_H_x[cind] = fj3[j]*_H_x[cind]
						  + fj2[j]*_CourantFactor*curle + fi1[i]*_I_Hx[cind];
  				
			}
		}
		// Calculate Hy
		for (auto j=0; j< _mesh->reg_num_nodes_y(); j++){
			for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);

				curle = _En_z[rind] - _En_z[cind];
				/*
				_I_Hy[cind] = _I_Hy[cind] + fj1[j]*curle;
				_H_y[cind] = fi3[i]*_H_y[cind]
						  + fi2[i]*_CourantFactor*(curle + _I_Hy[cind]);
				*/
				
				_I_Hy[cind] = _I_Hy[cind] + curle;
				_H_y[cind] = fi3[i]*_H_y[cind]
						  + fi2[i]*_CourantFactor*curle + fj1[j]*_I_Hy[cind];
				
			}
		}
		

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", _E_z[0]);
		//_simdata.add_data_at_index(n, "H_y", _H_y[0]);

		// update iteration count
		_current_iter++;

	}


}


void FDTDSimulation::run_3D(int num_iters){

	
	double pulse;
	double _srcx, _srcy, _srcz;

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else {
		end_iter = _current_iter + num_iters;
		if (end_iter > _num_iters) end_iter = _num_iters;
	}
	
	double curle, curlh, sourcemodval;
	unsigned int cind, lind, rind, uind, dind, chind, stind; // center, left, right, up, down, charm, strange
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;
		sourcemodval = _source_modulator.value(double(_current_iter)/_source_modulator_width*SIGNALGENERATOR_PI/2.0);
		//cout << endl;

		// update Dx field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){
				for (auto k=1; k<_mesh->reg_num_nodes_z()-1; k++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					//lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					//rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					//uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					//chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);
					

					curlh = (_H_z[cind] - _H_z[dind] - _H_y[cind] + _H_y[stind]);
					_I_Dnx[cind] = _I_Dnx[cind] + curlh;
					_Dn_x[cind] = gj3[j]*gk3[k]*_Dn_x[cind]
							  + gj2[j]*gk2[k]*_CourantFactor*(curlh + gi1[i]*_I_Dnx[cind])
							  - _CourantFactor*_dx*_current_density_x_fn(cind)*sourcemodval;
				}
			}
		}

		// update Dy field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){
				for (auto k=1; k<_mesh->reg_num_nodes_z()-1; k++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					//rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					//uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					//dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					//chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);
					

					curlh = (_H_x[cind] - _H_x[stind] - _H_z[cind] + _H_z[lind]);
					_I_Dny[cind] = _I_Dny[cind] + curlh;
					_Dn_y[cind] = gi3[i]*gk3[k]*_Dn_y[cind]
							  + gi2[i]*gk2[k]*_CourantFactor*(curlh + gj1[j]*_I_Dny[cind])
							  - _CourantFactor*_dx*_current_density_y_fn(cind)*sourcemodval;
				}
			}
		}

		// update Dz field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){
				for (auto k=1; k<_mesh->reg_num_nodes_z()-1; k++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					//rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					//uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					//chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					//stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);
					

					curlh = (_H_y[cind] - _H_y[lind] - _H_x[cind] + _H_x[dind]);
					_I_Dnz[cind] = _I_Dnz[cind] + curlh;
					_Dn_z[cind] = gi3[i]*gj3[j]*_Dn_z[cind]
							  + gi2[i]*gj2[j]*_CourantFactor*(curlh + gk1[k]*_I_Dnz[cind])
							  - _CourantFactor*_dx*_current_density_z_fn(cind)*sourcemodval;
				}
			}
		}

		//cout << "about to impose source" << endl;
		// impose sources
		for (auto i=0; i<_signals.size(); i++){
			pulse = _signals.at(i).value(_tcur)*_modulators.at(i).value(_tcur) * sourcemodval;
			_srcx = _signals.at(i).xloc();
			_srcy = _signals.at(i).yloc();
			_srcz = _signals.at(i).zloc();
			//for (auto j=-5; j<6; j++) _Dn_z[_mesh->nearest_node(_srcx, _srcy, _srcz+j*_dx)] = 0;
			_Dn_z[_mesh->nearest_node(_srcx, _srcy, _srcz)] = pulse;
			//cout << "Got here. res is " << _mesh->res() << endl;
			//cout << "Got here. source location: " << _srcx << ", " << _srcy << ", " << _srcz << " src value: " << pulse << endl;
		}
		//cout << "imposed source" << endl;
		
		// calculate Ez field
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				for (auto k=0; k<_mesh->reg_num_nodes_z(); k++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);

					// with pml, conductivity, and single pole
					_En_x[cind] = (_Dn_x[cind] - _I_Ex[cind] - exp(-_dt*_permittivity_single_pole_freq_fn(cind))*_S_nm1_x[cind])/
								  (_rel_permittivity_fn(cind) + _conductivity_fn(cind)*_dt/_eps0 + _permittivity_single_pole_numerator_fn(cind)*_dt);
					_En_y[cind] = (_Dn_y[cind] - _I_Ey[cind] - exp(-_dt*_permittivity_single_pole_freq_fn(cind))*_S_nm1_y[cind])/
								  (_rel_permittivity_fn(cind) + _conductivity_fn(cind)*_dt/_eps0 + _permittivity_single_pole_numerator_fn(cind)*_dt);
					_En_z[cind] = (_Dn_z[cind] - _I_Ez[cind] - exp(-_dt*_permittivity_single_pole_freq_fn(cind))*_S_nm1_z[cind])/
								  (_rel_permittivity_fn(cind) + _conductivity_fn(cind)*_dt/_eps0 + _permittivity_single_pole_numerator_fn(cind)*_dt);

					_I_Ex[cind] = _I_Ex[cind] + _conductivity_fn(cind)*_dt/_eps0*_En_x[cind];
					_I_Ey[cind] = _I_Ey[cind] + _conductivity_fn(cind)*_dt/_eps0*_En_y[cind];
					_I_Ez[cind] = _I_Ez[cind] + _conductivity_fn(cind)*_dt/_eps0*_En_z[cind];
					

					_S_n_x[cind] = exp(-_permittivity_single_pole_freq_fn(cind)*_dt)*_S_nm1_x[cind] + _permittivity_single_pole_numerator_fn(cind)*_dt*_En_x[cind];
					_S_n_y[cind] = exp(-_permittivity_single_pole_freq_fn(cind)*_dt)*_S_nm1_y[cind] + _permittivity_single_pole_numerator_fn(cind)*_dt*_En_y[cind];
					_S_n_z[cind] = exp(-_permittivity_single_pole_freq_fn(cind)*_dt)*_S_nm1_z[cind] + _permittivity_single_pole_numerator_fn(cind)*_dt*_En_z[cind];
					
					_S_nm1_x[cind] = _S_n_x[cind];
					_S_nm1_y[cind] = _S_n_y[cind];
					_S_nm1_z[cind] = _S_n_z[cind];

					_E_x[cind] = _En_x[cind] * sqrt(_mu0/_eps0);
					_E_y[cind] = _En_y[cind] * sqrt(_mu0/_eps0);
					_E_z[cind] = _En_z[cind] * sqrt(_mu0/_eps0); // this can be made faster
				}
			}
		}
		//cout << "data: " << _Dn_z[_mesh->nearest_node(_srcx, _srcy, _srcz)] << endl;

		
		// set E at boundaries to zero for pml
		// right and left boundaries
		for (auto i=0; i<_mesh->reg_num_nodes_z(); i++){
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				rind = _mesh->reg_inds_to_glob_ind(_mesh->reg_num_nodes_x()-1, j, i);
				lind = _mesh->reg_inds_to_glob_ind(0, j, i);
				_En_x[rind] = 0.0;
				_En_x[lind] = 0.0;
				_En_y[rind] = 0.0;
				_En_y[lind] = 0.0;
				_En_z[rind] = 0.0;
				_En_z[lind] = 0.0;
			}
		}
		// top and bottom boundaries
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			for (auto j=0; j<_mesh->reg_num_nodes_z(); j++){
				dind = _mesh->reg_inds_to_glob_ind(i, 0, j);
				uind = _mesh->reg_inds_to_glob_ind(i, _mesh->reg_num_nodes_y()-1, j);
				_En_x[dind] = 0.0;
				_En_x[uind] = 0.0;
				_En_y[dind] = 0.0;
				_En_y[uind] = 0.0;
				_En_z[dind] = 0.0;
				_En_z[uind] = 0.0;
			}
		}
		// charm and strange boundaries
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				chind = _mesh->reg_inds_to_glob_ind(i, j, _mesh->reg_num_nodes_z()-1);
				stind = _mesh->reg_inds_to_glob_ind(i, j, 0);
				_En_x[chind] = 0.0;
				_En_x[stind] = 0.0;
				_En_y[chind] = 0.0;
				_En_y[stind] = 0.0;
				_En_z[chind] = 0.0;
				_En_z[stind] = 0.0;
			}
		}
		
		
		

		// Calculate Hx
		for (auto k=0; k<_mesh->reg_num_nodes_z()-1; k++){
			for (auto j=0; j<_mesh->reg_num_nodes_y()-1; j++){
				for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					//lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					//rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					//dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					//stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);

					curle = _En_y[chind] - _En_y[cind] - _En_z[uind] + _En_z[cind];
					_I_Hx[cind] = _I_Hx[cind] + curle;
					_H_x[cind] = fj3[j]*fk3[k]*_H_x[cind]
							  + fj2[j]*fk2[k]*_CourantFactor*(curle + fi1[i]*_I_Hx[cind]);
				}
			}
		}
		// Calculate Hy
		for (auto k=0; k<_mesh->reg_num_nodes_z()-1; k++){
			for (auto j=0; j< _mesh->reg_num_nodes_y(); j++){
				for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					//lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					//uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					//dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					//stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);

					curle = _En_z[rind] - _En_z[cind] - _En_x[chind] + _En_x[cind];
					_I_Hy[cind] = _I_Hy[cind] + curle;
					_H_y[cind] = fi3[i]*fk3[k]*_H_y[cind]
							  + fi2[i]*fk2[k]*_CourantFactor*(curle + fj1[j]*_I_Hy[cind]);
				}
			}
		}
		// Calculate Hz
		for (auto k=0; k<_mesh->reg_num_nodes_z(); k++){
			for (auto j=0; j< _mesh->reg_num_nodes_y()-1; j++){
				for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
					cind = _mesh->reg_inds_to_glob_ind(i, j, k);
					//lind = _mesh->reg_inds_to_glob_ind(i-1, j, k);
					rind = _mesh->reg_inds_to_glob_ind(i+1, j, k);
					uind = _mesh->reg_inds_to_glob_ind(i, j+1, k);
					//dind = _mesh->reg_inds_to_glob_ind(i, j-1, k);
					//chind = _mesh->reg_inds_to_glob_ind(i, j, k+1);
					//stind = _mesh->reg_inds_to_glob_ind(i, j, k-1);


					curle = _En_x[uind] - _En_x[cind] - _En_y[rind] + _En_y[cind];
					_I_Hz[cind] = _I_Hz[cind] + curle;
					_H_z[cind] = fi3[i]*fj3[j]*_H_z[cind]
							  + fi2[i]*fj2[j]*_CourantFactor*(curle + fk1[k]*_I_Hz[cind]);
				}
			}
		}
		

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", _E_z[0]);
		//if (_E_z == nullptr) cout << "Ruh roh nullptr" << endl;
		cout << "near source : En_z:" << _En_z[_mesh->nearest_node(_srcx+_dx, _srcy, _srcz)] << "\tDn_z:" << _Dn_z[_mesh->nearest_node(_srcx+_dx, _srcy, _srcz)] << endl;

		//_simdata.add_data_at_index(n, "H_y", _H_y[0]);

		// update iteration count
		_current_iter++;

	}


}


void FDTDSimulation::preRunCheck(){
	
}

void FDTDSimulation::prepareModulator(){
	_source_modulator.set_tanh();
	_source_modulator_width = 20;
}

void FDTDSimulation::allocate_fields(){

	cout << "allocating fields" << endl;

	if (_mesh->num_dims() == 1){

		_Dn_z = new double[_mesh->nodecount()];
		_I_Hy = new double[_mesh->nodecount()];
		_En_z = new double[_mesh->nodecount()];
		_E_z = new double[_mesh->nodecount()];
		_H_y = new double[_mesh->nodecount()];
		_I_Ez = new double[_mesh->nodecount()];
		_S_n_z = new double[_mesh->nodecount()];
		_S_nm1_z = new double[_mesh->nodecount()];

		for (auto i=0; i<_mesh->nodecount(); i++){
			_Dn_z[i] = 0.0;
			_I_Hy[i] = 0.0;
			_H_y[i] = 0.0;
			_En_z[i] = 0.0;
			_E_z[i] = 0.0;
			_I_Ez[i] = 0.0;
			_S_n_z[i] = 0.0;
			_S_nm1_z[i] = 0.0;
		}

	}
	else if (_mesh->num_dims() == 2){

		_Dn_z = new double[_mesh->nodecount()];
		_I_Hx = new double[_mesh->nodecount()];
		_I_Hy = new double[_mesh->nodecount()];
		_En_z = new double[_mesh->nodecount()];
		_E_z = new double[_mesh->nodecount()];
		_H_y = new double[_mesh->nodecount()];
		_H_x = new double[_mesh->nodecount()];
		_I_Ez = new double[_mesh->nodecount()];
		_S_n_z = new double[_mesh->nodecount()];
		_S_nm1_z = new double[_mesh->nodecount()];

		for (auto i=0; i<_mesh->nodecount(); i++){
			_Dn_z[i] = 0.0;
			_I_Hx[i] = 0.0;
			_I_Hy[i] = 0.0;
			_H_x[i] = 0.0;
			_H_y[i] = 0.0;
			_En_z[i] = 0.0;
			_E_z[i] = 0.0;
			_I_Ez[i] = 0.0;
			_S_n_z[i] = 0.0;
			_S_nm1_z[i] = 0.0;
		}

	}
	else if (_mesh->num_dims() == 3){
		_Dn_x = new double[_mesh->nodecount()];
		_Dn_y = new double[_mesh->nodecount()];
		_Dn_z = new double[_mesh->nodecount()];
		_I_Dnx = new double[_mesh->nodecount()];
		_I_Dny = new double[_mesh->nodecount()];
		_I_Dnz = new double[_mesh->nodecount()];
		_I_Hx = new double[_mesh->nodecount()];
		_I_Hy = new double[_mesh->nodecount()];
		_I_Hz = new double[_mesh->nodecount()];
		_En_x = new double[_mesh->nodecount()];
		_En_y = new double[_mesh->nodecount()];
		_En_z = new double[_mesh->nodecount()];
		_E_x = new double[_mesh->nodecount()];
		_E_y = new double[_mesh->nodecount()];
		_E_z = new double[_mesh->nodecount()];
		_H_x = new double[_mesh->nodecount()];
		_H_y = new double[_mesh->nodecount()];
		_H_z = new double[_mesh->nodecount()];

		_I_Ex = new double[_mesh->nodecount()];
		_I_Ey = new double[_mesh->nodecount()];
		_I_Ez = new double[_mesh->nodecount()];
		_S_n_x = new double[_mesh->nodecount()];
		_S_n_y = new double[_mesh->nodecount()];
		_S_n_z = new double[_mesh->nodecount()];
		_S_nm1_x = new double[_mesh->nodecount()];
		_S_nm1_y = new double[_mesh->nodecount()];
		_S_nm1_z = new double[_mesh->nodecount()];



		for (auto i=0; i<_mesh->nodecount(); i++){
			_Dn_x[i] = 0.0;
			_Dn_y[i] = 0.0;
			_Dn_z[i] = 0.0;
			_I_Dnx[i] = 0.0;
			_I_Dny[i] = 0.0;
			_I_Dnz[i] = 0.0;
			_I_Hx[i] = 0.0;
			_I_Hy[i] = 0.0;
			_I_Hz[i] = 0.0;
			_H_x[i] = 0.0;
			_H_y[i] = 0.0;
			_H_z[i] = 0.0;
			_En_x[i] = 0.0;
			_En_y[i] = 0.0;
			_En_z[i] = 0.0;
			_E_x[i] = 0.0;
			_E_y[i] = 0.0;
			_E_z[i] = 0.0;

			_I_Ex[i] = 0.0;
			_I_Ey[i] = 0.0;
			_I_Ez[i] = 0.0;
			_S_n_x[i] = 0.0;
			_S_n_y[i] = 0.0;
			_S_n_z[i] = 0.0;
			_S_nm1_x[i] = 0.0;
			_S_nm1_y[i] = 0.0;
			_S_nm1_z[i] = 0.0;
		}

	}
	

}

void FDTDSimulation::allocate_coeffs(){

	
	// deal with relative permittivity
	if (_rel_permittivity_fn == nullptr){
		if (_rel_permittivity == nullptr) _rel_permittivity_fn = [](unsigned int i)->double{return 1.0;};
		else _rel_permittivity_fn = [this](unsigned int i)->double{return this->_rel_permittivity[i];};
	}
	
	// conductivity
	if (_conductivity_fn == nullptr){
		if (_conductivity == nullptr) _conductivity_fn = [](unsigned int i)->double{return 0.0;};
		else _conductivity_fn = [this](unsigned int i)->double{return this->_conductivity[i];};
	}

	// single pole material frequency (Hz)
	if (_permittivity_single_pole_freq_fn == nullptr){
		if (_permittivity_single_pole_freq == nullptr) _permittivity_single_pole_freq_fn = [](unsigned int i)->double{return 0.0;};
		else _permittivity_single_pole_freq_fn = [this](unsigned int i)->double{return this->_permittivity_single_pole_freq[i];};
	}

	// single pole material numerator
	if (_permittivity_single_pole_numerator_fn == nullptr){
		if (_permittivity_single_pole_numerator == nullptr) _permittivity_single_pole_numerator_fn = [](unsigned int i)->double{return 0.0;};
		else _permittivity_single_pole_numerator_fn = [this](unsigned int i)->double{return this->_permittivity_single_pole_numerator[i];};
	}

	// deal with relative permeability
	if (_rel_permeability_fn == nullptr){
		if (_rel_permeability == nullptr){

		}
	}

	// deal with current density
	if (_current_density_x_fn == nullptr){
		if (_current_density_x == nullptr) _current_density_x_fn = [](unsigned int i)->double{return 0.0;};
		else _current_density_x_fn = [this](unsigned int i)->double{return this->_current_density_x[i];};
	}

	if (_current_density_y_fn == nullptr){
		if (_current_density_y == nullptr) _current_density_y_fn = [](unsigned int i)->double{return 0.0;};
		else _current_density_y_fn = [this](unsigned int i)->double{return this->_current_density_y[i];};
	}

	if (_current_density_z_fn == nullptr){
		if (_current_density_z == nullptr) _current_density_z_fn = [](unsigned int i)->double{return 0.0;};
		else _current_density_z_fn = [this](unsigned int i)->double{return this->_current_density_z[i];};
	}

}

void FDTDSimulation::allocate_PML(){
	// PML related

	if (_mesh->num_dims() == 1){

		float xnum, xd, xxn, xn;
		gi1 = new double[_mesh->reg_num_nodes_x()];
		gi2 = new double[_mesh->reg_num_nodes_x()];
		gi3 = new double[_mesh->reg_num_nodes_x()];
		fi1 = new double[_mesh->reg_num_nodes_x()];
		fi2 = new double[_mesh->reg_num_nodes_x()];
		fi3 = new double[_mesh->reg_num_nodes_x()];
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			gi2[i] = 1.0;
			gi3[i] = 1.0;
			fi1[i] = 0.0;
			fi2[i] = 1.0;
			fi3[i] = 1.0;
		}


		if (_nPML == 0) return;

		// 
		for (auto i=0; i<= _nPML; i++){
			xnum = _nPML-i;
			xd = _nPML;

			// left and right boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gi1[i] = xn;
			gi2[i] = 1.0/(1.0+xn);
			gi3[i] = (1.0-xn)/(1.0+xn);
			gi1[_mesh->reg_num_nodes_x()-1-i] = xn;
			gi2[_mesh->reg_num_nodes_x()-1-i] = 1.0/(1.0+xn);
			gi3[_mesh->reg_num_nodes_x()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fi1[i] = xn;
			fi2[i] = 1.0/(1.0+xn);
			fi3[i] = (1.0-xn)/(1.0+xn);
			fi1[_mesh->reg_num_nodes_x()-2-i] = xn;
			fi2[_mesh->reg_num_nodes_x()-2-i] = 1.0/(1.0+xn);
			fi3[_mesh->reg_num_nodes_x()-2-i] = (1.0-xn)/(1.0+xn);
		}
	}
	else if(_mesh->num_dims() == 2){

		float xnum, xd, xxn, xn;
		gi1 = new double[_mesh->reg_num_nodes_x()];
		gi2 = new double[_mesh->reg_num_nodes_x()];
		gi3 = new double[_mesh->reg_num_nodes_x()];
		gj1 = new double[_mesh->reg_num_nodes_y()];
		gj2 = new double[_mesh->reg_num_nodes_y()];
		gj3 = new double[_mesh->reg_num_nodes_y()];
		fi1 = new double[_mesh->reg_num_nodes_x()];
		fi2 = new double[_mesh->reg_num_nodes_x()];
		fi3 = new double[_mesh->reg_num_nodes_x()];
		fj1 = new double[_mesh->reg_num_nodes_y()];
		fj2 = new double[_mesh->reg_num_nodes_y()];
		fj3 = new double[_mesh->reg_num_nodes_y()];
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			gi1[i] = 1.0;
			gi2[i] = 1.0;
			gi3[i] = 1.0;
			fi1[i] = 0.0;
			fi2[i] = 1.0;
			fi3[i] = 1.0;
		}
		for (auto i=0; i<_mesh->reg_num_nodes_y(); i++){
			gj1[i] = 1.0;
			gj2[i] = 1.0;
			gj3[i] = 1.0;
			fj1[i] = 0.0;
			fj2[i] = 1.0;
			fj3[i] = 1.0;
		}


		if (_nPML == 0) return;

		// 
		for (auto i=0; i<= _nPML; i++){
			xnum = _nPML-i;
			xd = _nPML;

			// left and right boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gi1[i] = xn;
			gi2[i] = 1.0/(1.0+xn);
			gi3[i] = (1.0-xn)/(1.0+xn);
			gi1[_mesh->reg_num_nodes_x()-1-i] = xn;
			gi2[_mesh->reg_num_nodes_x()-1-i] = 1.0/(1.0+xn);
			gi3[_mesh->reg_num_nodes_x()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fi1[i] = xn;
			fi2[i] = 1.0/(1.0+xn);
			fi3[i] = (1.0-xn)/(1.0+xn);
			fi1[_mesh->reg_num_nodes_x()-2-i] = xn;
			fi2[_mesh->reg_num_nodes_x()-2-i] = 1.0/(1.0+xn);
			fi3[_mesh->reg_num_nodes_x()-2-i] = (1.0-xn)/(1.0+xn);


			// top and bottom boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gj1[i] = xn;
			gj2[i] = 1.0/(1.0+xn);
			gj3[i] = (1.0-xn)/(1.0+xn);
			gj1[_mesh->reg_num_nodes_y()-1-i] = xn;
			gj2[_mesh->reg_num_nodes_y()-1-i] = 1.0/(1.0+xn);
			gj3[_mesh->reg_num_nodes_y()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fj1[i] = xn;
			fj2[i] = 1.0/(1.0+xn);
			fj3[i] = (1.0-xn)/(1.0+xn);
			fj1[_mesh->reg_num_nodes_y()-2-i] = xn;
			fj2[_mesh->reg_num_nodes_y()-2-i] = 1.0/(1.0+xn);
			fj3[_mesh->reg_num_nodes_y()-2-i] = (1.0-xn)/(1.0+xn);
		}

	}
	else if(_mesh->num_dims() == 3){

		float xnum, xd, xxn, xn;
		gi1 = new double[_mesh->reg_num_nodes_x()];
		gi2 = new double[_mesh->reg_num_nodes_x()];
		gi3 = new double[_mesh->reg_num_nodes_x()];
		gj1 = new double[_mesh->reg_num_nodes_y()];
		gj2 = new double[_mesh->reg_num_nodes_y()];
		gj3 = new double[_mesh->reg_num_nodes_y()];
		gk1 = new double[_mesh->reg_num_nodes_z()];
		gk2 = new double[_mesh->reg_num_nodes_z()];
		gk3 = new double[_mesh->reg_num_nodes_z()];
		fi1 = new double[_mesh->reg_num_nodes_x()];
		fi2 = new double[_mesh->reg_num_nodes_x()];
		fi3 = new double[_mesh->reg_num_nodes_x()];
		fj1 = new double[_mesh->reg_num_nodes_y()];
		fj2 = new double[_mesh->reg_num_nodes_y()];
		fj3 = new double[_mesh->reg_num_nodes_y()];
		fk1 = new double[_mesh->reg_num_nodes_z()];
		fk2 = new double[_mesh->reg_num_nodes_z()];
		fk3 = new double[_mesh->reg_num_nodes_z()];
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			gi1[i] = 1.0;
			gi2[i] = 1.0;
			gi3[i] = 1.0;
			fi1[i] = 0.0;
			fi2[i] = 1.0;
			fi3[i] = 1.0;
		}
		for (auto i=0; i<_mesh->reg_num_nodes_y(); i++){
			gj1[i] = 1.0;
			gj2[i] = 1.0;
			gj3[i] = 1.0;
			fj1[i] = 0.0;
			fj2[i] = 1.0;
			fj3[i] = 1.0;
		}
		for (auto i=0; i<_mesh->reg_num_nodes_z(); i++){
			gk1[i] = 1.0;
			gk2[i] = 1.0;
			gk3[i] = 1.0;
			fk1[i] = 0.0;
			fk2[i] = 1.0;
			fk3[i] = 1.0;
		}


		if (_nPML == 0) return;

		// 
		for (auto i=0; i<= _nPML; i++){
			xnum = _nPML-i;
			xd = _nPML;

			// left and right boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gi1[i] = xn;
			gi2[i] = 1.0/(1.0+xn);
			gi3[i] = (1.0-xn)/(1.0+xn);
			gi1[_mesh->reg_num_nodes_x()-1-i] = xn;
			gi2[_mesh->reg_num_nodes_x()-1-i] = 1.0/(1.0+xn);
			gi3[_mesh->reg_num_nodes_x()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fi1[i] = xn;
			fi2[i] = 1.0/(1.0+xn);
			fi3[i] = (1.0-xn)/(1.0+xn);
			fi1[_mesh->reg_num_nodes_x()-2-i] = xn;
			fi2[_mesh->reg_num_nodes_x()-2-i] = 1.0/(1.0+xn);
			fi3[_mesh->reg_num_nodes_x()-2-i] = (1.0-xn)/(1.0+xn);


			// top and bottom boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gj1[i] = xn;
			gj2[i] = 1.0/(1.0+xn);
			gj3[i] = (1.0-xn)/(1.0+xn);
			gj1[_mesh->reg_num_nodes_y()-1-i] = xn;
			gj2[_mesh->reg_num_nodes_y()-1-i] = 1.0/(1.0+xn);
			gj3[_mesh->reg_num_nodes_y()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fj1[i] = xn;
			fj2[i] = 1.0/(1.0+xn);
			fj3[i] = (1.0-xn)/(1.0+xn);
			fj1[_mesh->reg_num_nodes_y()-2-i] = xn;
			fj2[_mesh->reg_num_nodes_y()-2-i] = 1.0/(1.0+xn);
			fj3[_mesh->reg_num_nodes_y()-2-i] = (1.0-xn)/(1.0+xn);
		
			// charm and strange boundaries
			xxn = xnum/xd;
			xn = 0.33*xxn*xxn*xxn;
			gk1[i] = xn;
			gk2[i] = 1.0/(1.0+xn);
			gk3[i] = (1.0-xn)/(1.0+xn);
			gk1[_mesh->reg_num_nodes_z()-1-i] = xn;
			gk2[_mesh->reg_num_nodes_z()-1-i] = 1.0/(1.0+xn);
			gk3[_mesh->reg_num_nodes_z()-1-i] = (1.0-xn)/(1.0+xn);

			xxn = (xnum-0.5)/xd;
			xn = 0.25*xxn*xxn*xxn;
			fk1[i] = xn;
			fk2[i] = 1.0/(1.0+xn);
			fk3[i] = (1.0-xn)/(1.0+xn);
			fk1[_mesh->reg_num_nodes_z()-2-i] = xn;
			fk2[_mesh->reg_num_nodes_z()-2-i] = 1.0/(1.0+xn);
			fk3[_mesh->reg_num_nodes_z()-2-i] = (1.0-xn)/(1.0+xn);

		}
	}
	else{
		cout << "WHY YOU NO DEFINE NUMBER OF DIMENSIONS IN MESH?" << endl;
		throw -1;
	}
	
}

void FDTDSimulation::allocate_simdata(){

	// set up the simulation data holder
	_simdata.bind_mesh(*_mesh);
	_simdata.add_field("E_z");
	//_simdata.add_field("H_y");
	_simdata.set_time_span(0.0, _dt, _num_iters*_dt);
	_simdata.print_summary();
}

