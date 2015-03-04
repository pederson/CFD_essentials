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

	
	double * D_z = new double[_mesh.nodecount()];
	double * I_Hx = new double[_mesh.nodecount()];
	double * I_Hy = new double[_mesh.nodecount()];
	double * E_z = new double[_mesh.nodecount()];
	double * H_y = new double[_mesh.nodecount()];
	double * H_x = new double[_mesh.nodecount()];
	double * gaz = new double[_mesh.nodecount()];
	//double * gbz = new double[_mesh.nodecount()];
	for (auto i=0; i<_mesh.nodecount(); i++){
		D_z[i] = 0.0;
		I_z[i] = 0.0;
		I_Hx[i] = 0.0;
		I_Hy[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
		gaz[i] = 1.0/_rel_permittivity[i];		// gaz = 1/(epsilon+(sigma*dt/eps0));
		//gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

	// PML related
	float npml, xnum, xd, xxn, xn;
	npml = 8;
	double * gi2 = new double[_mesh.reg_num_nodes_x()];
	double * gi3 = new double[_mesh.reg_num_nodes_x()];
	double * gj2 = new double[_mesh.reg_num_nodes_y()];
	double * gj3 = new double[_mesh.reg_num_nodes_y()];
	double * fi1 = new double[_mesh.reg_num_nodes_x()];
	double * fi2 = new double[_mesh.reg_num_nodes_x()];
	double * fi3 = new double[_mesh.reg_num_nodes_x()];
	double * fj1 = new double[_mesh.reg_num_nodes_y()];
	double * fj2 = new double[_mesh.reg_num_nodes_y()];
	double * fj3 = new double[_mesh.reg_num_nodes_y()];
	for (auto i=0; i<_mesh.reg_num_nodes_x(); i++){
		gi2[i] = 1.0;
		gi3[i] = 1.0;
		fi1[i] = 0.0;
		fi2[i] = 1.0;
		fi3[i] = 1.0;
	}
	for (auto i=0; i<_mesh.reg_num_nodes_y(); i++){
		gj2[i] = 1.0;
		gj3[i] = 1.0;
		fj1[i] = 0.0;
		fj2[i] = 1.0;
		fj3[i] = 1.0;
	}
	// 
	for (auto i=0; i<= npml; i++){
		xnum = npml-i;
		xd = npml;

		// left and right boundaries
		xxn = xnum/xd;
		xn = 0.33*xxn*xxn*xxn;
		gi2[i] = 1.0/(1.0+xn);
		gi3[i] = (1.0-xn)/(1.0+xn);
		gi2[_mesh.reg_num_nodes_x()-1-i] = 1.0/(1.0+xn);
		gi3[_mesh.reg_num_nodes_x()-1-i] = (1.0-xn)/(1.0+xn);

		xxn = (xnum-0.5)/xd;
		xn = 0.25*xxn*xxn*xxn;
		fi1[i] = xn;
		fi2[i] = 1.0/(1.0+xn);
		fi3[i] = (1.0-xn)/(1.0+xn);
		fi1[_mesh.reg_num_nodes_x()-2-i] = xn;
		fi2[_mesh.reg_num_nodes_x()-2-i] = 1.0/(1.0+xn);
		fi3[_mesh.reg_num_nodes_x()-2-i] = (1.0-xn)/(1.0+xn);


		// top and bottom boundaries
		xxn = xnum/xd;
		xn = 0.33*xxn*xxn*xxn;
		gj2[i] = 1.0/(1.0+xn);
		gj3[i] = (1.0-xn)/(1.0+xn);
		gj2[_mesh.reg_num_nodes_y()-1-i] = 1.0/(1.0+xn);
		gj3[_mesh.reg_num_nodes_y()-1-i] = (1.0-xn)/(1.0+xn);

		xxn = (xnum-0.5)/xd;
		xn = 0.25*xxn*xxn*xxn;
		fj1[i] = xn;
		fj2[i] = 1.0/(1.0+xn);
		fj3[i] = (1.0-xn)/(1.0+xn);
		fj1[_mesh.reg_num_nodes_y()-2-i] = xn;
		fj2[_mesh.reg_num_nodes_y()-2-i] = 1.0/(1.0+xn);
		fj3[_mesh.reg_num_nodes_y()-2-i] = (1.0-xn)/(1.0+xn);


	}

	double CourantFactor = 0.5;
	double tcur, pulse, t0=10.0, spread = 12.0, srcfreq=1.0e+9, curle;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=0; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;


		// update D field
		for (auto i=1; i<_mesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh.reg_num_nodes_y()-1; j++){

				cind = _mesh.reg_inds_to_glob_ind(i, j);
				lind = _mesh.reg_inds_to_glob_ind(i-1, j);
				rind = _mesh.reg_inds_to_glob_ind(i+1, j);
				uind = _mesh.reg_inds_to_glob_ind(i, j+1);
				dind = _mesh.reg_inds_to_glob_ind(i, j-1);

				
				D_z[cind] = gi3[i]*gj3[j]*D_z[cind]
						  + gi2[i]*gj2[j]*CourantFactor*(H_y[cind] - H_y[lind] - H_x[cind] + H_x[dind]);

				
			}
		}


		// impose a gaussian source (hard-coded)
		//pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		pulse = sin(2*VX_PI*srcfreq*n*dt); // oscillatory source
		D_z[_mesh.nearest_node(0.5, 0.5)] = pulse;

		
		// calculate Ez field
		for (auto i=0; i<_mesh.reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<_mesh.reg_num_nodes_y(); j++){
				cind = _mesh.reg_inds_to_glob_ind(i, j);

				// with pml
				E_z[cind] = gaz[cind] * D_z[cind];
				
			}
		}
		// set Ez at edges to zero for pml
		for (auto j=0; j<_mesh.reg_num_nodes_y(); j++){
			lind = _mesh.reg_inds_to_glob_ind(0, j);
			rind = _mesh.reg_inds_to_glob_ind(_mesh.reg_num_nodes_x()-1, j);
			E_z[lind] = 0.0;
			E_z[rind] = 0.0;
		}
		for (auto i=0; i<_mesh.reg_num_nodes_x(); i++){
			uind = _mesh.reg_inds_to_glob_ind(i, _mesh.reg_num_nodes_y()-1);
			dind = _mesh.reg_inds_to_glob_ind(i, 0);
			E_z[uind] = 0.0;
			E_z[dind] = 0.0;
		}
		

		// Calculate Hx
		for (auto j=0; j<_mesh.reg_num_nodes_y()-1; j++){
			for (auto i=0; i<_mesh.reg_num_nodes_x(); i++){
				cind = _mesh.reg_inds_to_glob_ind(i, j);
				uind = _mesh.reg_inds_to_glob_ind(i, j+1);

				curle = E_z[cind] - E_z[uind];
				I_Hx[cind] = I_Hx[cind] + fi1[i]*curle;
				H_x[cind] = fj3[j]*H_x[cind]
						  + fj2[j]*CourantFactor*(curle + I_Hx[cind]);

			}
		}
		// Calculate Hy
		for (auto j=0; j< _mesh.reg_num_nodes_y(); j++){
			for (auto i=0; i<_mesh.reg_num_nodes_x()-1; i++){
				cind = _mesh.reg_inds_to_glob_ind(i, j);
				rind = _mesh.reg_inds_to_glob_ind(i+1, j);

				curle = E_z[rind] - E_z[cind];
				I_Hy[cind] = I_Hy[cind] + fj1[j]*curle;
				H_y[cind] = fi3[i]*H_y[cind]
						  + fi2[i]*CourantFactor*(curle + I_Hy[cind]);

			}
		}
		

		// fill in the simdata for this time step
		simdata.add_data_at_index(n, "E_z", E_z[0]);
		simdata.add_data_at_index(n, "H_y", H_y[0]);

		tcur += dt;
	}

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
	delete[] I_Hx;
	delete[] I_Hy;
	delete[] E_z;
	delete[] H_y;
	delete[] H_x;
	delete[] gaz;

	delete[] gi2;
	delete[] gi3;
	delete[] gj2;
	delete[] gj3;
	delete[] fi1;
	delete[] fi2;
	delete[] fi3;
	delete[] fj1;
	delete[] fj2;
	delete[] fj3;

}

