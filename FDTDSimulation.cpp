#include "FDTDSimulation.hpp"


// mutators
void FDTDSimulation::set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=0);

void FDTDSimulation::add_gaussian_source(double t0, double spread, double xloc, double yloc=0, double zloc=0);

void FDTDSimulation::add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc=0, double zloc=0);

void FDTDSimulation::bind_rel_permittivity(const double * rel_permittivity){
	_rel_permittivity = rel_permittivity;
}

void FDTDSimulation::set_view_results(bool opt){
	_visualize_results = opt;
}

void FDTDSimulation::set_output_HDF5(bool opt, std::string outname){
	_output_HDF5 = opt;
	if (outname.compare("") == 0){
		_output_HDF5_name = "FDTD_output.h5";
	}
}

void FDTDSimulation::run(){

	// set up the simulation data holder
	_simdata.bind_mesh(_mesh);
	_simdata.add_field("E_z");
	_simdata.set_time_span(_tstart, _dt, _tstart);
	_simdata.print_summary();

	// calculate E using centered difference except for on the boundaries
	double * D_z = new double[_mesh.nodecount()];
	double * I_z = new double[_mesh.nodecount()];
	double * E_z = new double[_mesh.nodecount()];
	double * H_y = new double[_mesh.nodecount()];
	double * H_x = new double[_mesh.nodecount()];
	double * gaz = new double[_mesh.nodecount()];
	double * gbz = new double[_mesh.nodecount()];
	for (auto i=0; i<_mesh.nodecount(); i++){
		D_z[i] = 0.0;
		I_z[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
		gaz[i] = 1.0/_rel_permittivity[i];		// gaz = 1/(epsilon+(sigma*dt/eps0));
		gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

	double tcur;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=0; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;

		// update E field
		for (auto i=1; i<_mesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh.reg_num_nodes_y()-1; j++){
				cind = _mesh.reg_inds_to_glob_ind(i, j);
				lind = _mesh.reg_inds_to_glob_ind(i, j-1);
				rind = _mesh.reg_inds_to_glob_ind(i, j+1);
				uind = _mesh.reg_inds_to_glob_ind(i+1, j);
				dind = _mesh.reg_inds_to_glob_ind(i-1, j);

				//H_y[cind] = H_y[cind] + 0.5*(E_x[cind] - E_x[rind]);

				D_z[cind] = D_z[cind] 
						  + 0.5 * (H_y[cind] - H_y[lind] - H_x[cind] + H_x[dind]);
				E_z[cind] = gaz[cind] *(D_z[cind] - I_z[cind]);
				I_z[cind] = I_z[cind] + gbz[cind] * E_z[cind];
			}
		}



		// impose a gaussian source (hard-coded)
		//pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		pulse = sin(2*VX_PI*srcfreq*n*dt); // oscillatory source
		//E_z[paramesh.nodecount()/4] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.9)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.8)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.7)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.6)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.5)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.4)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.3)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.2)] = pulse;
		E_z[_mesh.nearest_node(0.25, 0.1)] = pulse;

		
		// update H field
		//cout << "calculating new Electric field" << endl;
		for (auto i=1; i<_mesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh.reg_num_nodes_y()-1; j++){
				cind = _mesh.reg_inds_to_glob_ind(i, j);
				lind = _mesh.reg_inds_to_glob_ind(i, j-1);
				rind = _mesh.reg_inds_to_glob_ind(i, j+1);
				uind = _mesh.reg_inds_to_glob_ind(i+1, j);
				dind = _mesh.reg_inds_to_glob_ind(i-1, j);

				//E_x[cind] += 0.5*(H_y[lind] - H_y[cind]);
				H_x[cind] = H_x[cind] - 0.5 * (E_z[uind] - E_z[cind]);
				H_y[cind] = H_y[cind] + 0.5 * (E_z[rind] - E_z[cind]);
			}
		}

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", E_z[0]);

		tcur += dt;
	}
	//*/

	if (_output_HDF5) _simdata.write_HDF5(_output_HDF5_name);

	if (_visualize_results){
		// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_colorby_field("E_z");
		simvis.set_color_interpolation(false);
		simvis.set_frequency_Hz(30);
		simvis.run();
		//*/
	}

	delete[] D_z;
	delete[] I_z;
	delete[] E_z;
	delete[] H_y;
	delete[] H_x;
	delete[] gaz;
	delete[] gbz;
}

