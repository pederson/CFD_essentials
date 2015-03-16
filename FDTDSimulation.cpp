#include "FDTDSimulation.hpp"

using namespace std;

FDTDSimulation::FDTDSimulation(){
	_Dn_z = nullptr;
	_I_Hx = nullptr;
	_I_Hy = nullptr;
	_En_z = nullptr;
	_E_z = nullptr;
	_H_y = nullptr;
	_H_x = nullptr;
	_gaz = nullptr;
	_jnx = nullptr; 
	_jny = nullptr;
	_jnz = nullptr;

	gi2 = nullptr;
	gi3 = nullptr;
	gj2 = nullptr;
	gj3 = nullptr;
	fi1 = nullptr;
	fi2 = nullptr;
	fi3 = nullptr;
	fj1 = nullptr;
	fj2 = nullptr;
	fj3 = nullptr;

	_rel_permittivity = nullptr;
	_rel_permeability = nullptr;
	_current_density_x = nullptr;
	_field_multiplier_x = nullptr;
	_current_density_y = nullptr;
	_field_multiplier_y = nullptr;
	_current_density_z = nullptr;
	_field_multiplier_z = nullptr;

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
	if (_Dn_z != nullptr) delete[] _Dn_z;
	if (_I_Hx != nullptr) delete[] _I_Hx;
	if (_I_Hy != nullptr) delete[] _I_Hy;
	if (_En_z != nullptr) delete[] _En_z;
	if (_E_z != nullptr) delete[] _E_z;
	if (_H_y != nullptr) delete[] _H_y;
	if (_H_x != nullptr) delete[] _H_x;
	if (_gaz != nullptr) delete[] _gaz;
	if (_jnx != nullptr) delete[] _jnx;
	if (_jny != nullptr) delete[] _jny;
	if (_jnz != nullptr) delete[] _jnz;

	if (gi2 != nullptr) delete[] gi2;
	if (gi3 != nullptr) delete[] gi3;
	if (gj2 != nullptr) delete[] gj2;
	if (gj3 != nullptr) delete[] gj3;
	if (fi1 != nullptr) delete[] fi1;
	if (fi2 != nullptr) delete[] fi2;
	if (fi3 != nullptr) delete[] fi3;
	if (fj1 != nullptr) delete[] fj1;
	if (fj2 != nullptr) delete[] fj2;
	if (fj3 != nullptr) delete[] fj3;
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
}

void FDTDSimulation::bind_rel_permittivity(const double * rel_permittivity){
	_rel_permittivity = rel_permittivity;
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



void FDTDSimulation::view_results(){

		// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_color_ramp(CRamp::MATLAB_PARULA);
		simvis.set_colorby_field("E_z");
		//simvis.set_colorby_field("H_y");
		simvis.set_color_alpha(_rel_permittivity);
		simvis.set_color_interpolation(false);
		simvis.set_frequency_Hz(30);
		simvis.run();

}

void FDTDSimulation::output_HDF5(std::string outname){
	_simdata.write_HDF5(outname);
}

void FDTDSimulation::run(int num_iters){
	_dt = _CourantFactor*_dx/_c0;


	if (!_is_allocated) {
		preRunCheck();
		allocate_fields();
		allocate_PML();
		allocate_simdata();
		_is_allocated = true;
	}
	//if (_is_allocated && _time_variable_coefficients) recalcCoeffs();

	cout << "dx: " << _dx << endl;
	cout << "Courant Factor: " << _CourantFactor << endl;
	cout << "dt: " << _dt << endl;


	if (_mesh->num_dims() == 1){
		cout << "1D FDTD is not yet plugged in!" << endl;
	}
	else if (_mesh->num_dims() == 2){
		run_2D(num_iters);
	}
	else if (_mesh->num_dims() == 3){
		cout << "3D FDTD is not yet implemented!" << endl;
	}

}

void FDTDSimulation::run_2D(int num_iters){

	cout << "starting the run 2D" << endl;
	double pulse;//, t0=10.0, spread = 12.0, srcfreq=3.0e+14/6.0;
	double _srcx, _srcy;// = 0.5e-6, _srcy = 3.0e-6;

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
						  - _CourantFactor*_dx*_jnz[cind]*sourcemodval;

				
			}
		}

		//cout << "about to impose source" << endl;
		// impose sources
		//pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		//pulse = sin(2*VX_PI*srcfreq*_tcur); // oscillatory source
		//D_z[_mesh->nearest_node(_srcx, _srcy)] = pulse;
		for (auto i=0; i<_signals.size(); i++){
			pulse = _signals.at(i).value(_tcur)*_modulators.at(i).value(_tcur) * sourcemodval;
			_srcx = _signals.at(i).xloc();
			_srcy = _signals.at(i).yloc();
			_Dn_z[_mesh->nearest_node(_srcx, _srcy)] = pulse;
		}
		//cout << "impose source" << endl;
		
		// calculate Ez field
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);

				// with pml
				_En_z[cind] = _gaz[cind] * _Dn_z[cind];
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
				_I_Hx[cind] = _I_Hx[cind] + fi1[i]*curle;
				_H_x[cind] = fj3[j]*_H_x[cind]
						  + fj2[j]*_CourantFactor*(curle + _I_Hx[cind]);

			}
		}
		// Calculate Hy
		for (auto j=0; j< _mesh->reg_num_nodes_y(); j++){
			for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);

				curle = _En_z[rind] - _En_z[cind];
				_I_Hy[cind] = _I_Hy[cind] + fj1[j]*curle;
				_H_y[cind] = fi3[i]*_H_y[cind]
						  + fi2[i]*_CourantFactor*(curle + _I_Hy[cind]);

			}
		}
		

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", _E_z[0]);
		_simdata.add_data_at_index(n, "H_y", _H_y[0]);

		// update iteration count
		_current_iter++;

	}


}


void FDTDSimulation::preRunCheck(){

	/*
	if (_current_density_x == nullptr && _current_density_y == nullptr && _current_density_z == nullptr){
		_include_current_density = false;
	}
	else {
		_include_current_density = true;		
	}
	*/
	
}

void FDTDSimulation::prepareModulator(){
	_source_modulator.set_tanh();
	_source_modulator_width = 20;
}

void FDTDSimulation::allocate_fields(){

	_Dn_z = new double[_mesh->nodecount()];
	_I_Hx = new double[_mesh->nodecount()];
	_I_Hy = new double[_mesh->nodecount()];
	_En_z = new double[_mesh->nodecount()];
	_E_z = new double[_mesh->nodecount()];
	_H_y = new double[_mesh->nodecount()];
	_H_x = new double[_mesh->nodecount()];
	_gaz = new double[_mesh->nodecount()];
	_jnx = new double[_mesh->nodecount()];
	_jny = new double[_mesh->nodecount()];
	_jnz = new double[_mesh->nodecount()];
	//gbz = new double[_mesh.nodecount()];


	for (auto i=0; i<_mesh->nodecount(); i++){
		_Dn_z[i] = 0.0;
		_I_Hx[i] = 0.0;
		_I_Hy[i] = 0.0;
		_H_x[i] = 0.0;
		_H_y[i] = 0.0;
		_En_z[i] = 0.0;
		_E_z[i] = 0.0;
		//gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

	// deal with relative permittivity
	if (_rel_permittivity == nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++) _gaz[i] = 1.0;
	}
	else{
		for (auto i=0; i<_mesh->nodecount(); i++) _gaz[i] = 1.0/_rel_permittivity[i];
	}

	// deal with relative permeability
	if (_rel_permeability == nullptr){

	}

	// deal with current density
	if (_current_density_x == nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++) _jnx[i] = 0.0;
	}
	else{
		for (auto i=0; i<_mesh->nodecount(); i++) _jnx[i] = _current_density_x[i];	
	}

	if (_current_density_y == nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++) _jny[i] = 0.0;
	}
	else{
		for (auto i=0; i<_mesh->nodecount(); i++) _jny[i] = _current_density_y[i];
	}

	if (_current_density_z == nullptr){
		for (auto i=0; i<_mesh->nodecount(); i++) _jnz[i] = 0.0;
	}
	else{
		for (auto i=0; i<_mesh->nodecount(); i++) _jnz[i] = _current_density_z[i];	
	}




}

void FDTDSimulation::allocate_PML(){

	// PML related
	float xnum, xd, xxn, xn;
	gi2 = new double[_mesh->reg_num_nodes_x()];
	gi3 = new double[_mesh->reg_num_nodes_x()];
	gj2 = new double[_mesh->reg_num_nodes_y()];
	gj3 = new double[_mesh->reg_num_nodes_y()];
	fi1 = new double[_mesh->reg_num_nodes_x()];
	fi2 = new double[_mesh->reg_num_nodes_x()];
	fi3 = new double[_mesh->reg_num_nodes_x()];
	fj1 = new double[_mesh->reg_num_nodes_y()];
	fj2 = new double[_mesh->reg_num_nodes_y()];
	fj3 = new double[_mesh->reg_num_nodes_y()];
	for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
		gi2[i] = 1.0;
		gi3[i] = 1.0;
		fi1[i] = 0.0;
		fi2[i] = 1.0;
		fi3[i] = 1.0;
	}
	for (auto i=0; i<_mesh->reg_num_nodes_y(); i++){
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
		gi2[i] = 1.0/(1.0+xn);
		gi3[i] = (1.0-xn)/(1.0+xn);
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
		gj2[i] = 1.0/(1.0+xn);
		gj3[i] = (1.0-xn)/(1.0+xn);
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

void FDTDSimulation::allocate_simdata(){

	// set up the simulation data holder
	_simdata.bind_mesh(*_mesh);
	_simdata.add_field("E_z");
	_simdata.add_field("H_y");
	_simdata.set_time_span(0.0, _dt, _num_iters*_dt);
	_simdata.print_summary();
}

