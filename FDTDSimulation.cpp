#include "FDTDSimulation.hpp"

using namespace std;

FDTDSimulation::FDTDSimulation(){
	D_z = nullptr;
	I_Hx = nullptr;
	I_Hy = nullptr;
	E_z = nullptr;
	H_y = nullptr;
	H_x = nullptr;
	gaz = nullptr;

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

	_dx = 0.0;
	_dt = 0.0;
	_CourantFactor = 0.5;
	_nPML = 8;
	_num_iters = 0;
	_current_iter = 0;
	_tcur = 0.0;


}

FDTDSimulation::~FDTDSimulation(){
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
void FDTDSimulation::set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers){

}

void FDTDSimulation::add_gaussian_source(double t0, double spread, double xloc, double yloc, double zloc){

}

void FDTDSimulation::add_sinusoidal_source(double freq_Hz, double phase, double xloc, double yloc, double zloc){

}

void FDTDSimulation::bind_mesh(const RegularMesh & mesh) {
	_mesh = &mesh;
	_dx = _mesh->res();
}

void FDTDSimulation::bind_rel_permittivity(const double * rel_permittivity){
	_rel_permittivity = rel_permittivity;
}

void FDTDSimulation::view_results(){

		// visualize the simulation
		simulation_visualixer simvis;
		simvis.bind_simulation(_simdata);
		simvis.set_colorby_field("E_z");
		simvis.set_color_alpha(_rel_permittivity);
		simvis.set_color_interpolation(false);
		simvis.set_frequency_Hz(30);
		simvis.run();

}

void FDTDSimulation::output_HDF5(std::string outname){
	if (outname.compare("") == 0){
		_output_HDF5_name = "FDTD_output.h5";
	}
	else _output_HDF5_name = outname;

	_simdata.write_HDF5(_output_HDF5_name);
}

void FDTDSimulation::run(int num_iters){
	_dt = _CourantFactor*_dx/_c0;

	preRunCheck();
	allocate_fields();
	allocate_PML();
	allocate_simdata();

	cout << "dx: " << _dx << endl;
	cout << "Courant Factor: " << _CourantFactor << endl;


	

	run_2D(num_iters);

}

void FDTDSimulation::run_2D(int num_iters){

	cout << "starting the run" << endl;

	unsigned int end_iter;
	if (num_iters==-1) end_iter = _num_iters;
	else end_iter = _current_iter + num_iters;
	
	double pulse, t0=10.0, spread = 12.0, srcfreq=1.0e+9, curle;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=_current_iter; n<end_iter; n++){
		_tcur += _dt;

		cout << "on time step " << _current_iter << "/" << _num_iters-1 << "\r" << flush;


		// update D field
		for (auto i=1; i<_mesh->reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<_mesh->reg_num_nodes_y()-1; j++){

				cind = _mesh->reg_inds_to_glob_ind(i, j);
				lind = _mesh->reg_inds_to_glob_ind(i-1, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);
				uind = _mesh->reg_inds_to_glob_ind(i, j+1);
				dind = _mesh->reg_inds_to_glob_ind(i, j-1);

				
				D_z[cind] = gi3[i]*gj3[j]*D_z[cind]
						  + gi2[i]*gj2[j]*_CourantFactor*(H_y[cind] - H_y[lind] - H_x[cind] + H_x[dind]);

				
			}
		}


		// impose sources
		//pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		pulse = sin(2*VX_PI*srcfreq*n*_dt); // oscillatory source
		D_z[_mesh->nearest_node(0.5, 0.5)] = pulse;

		
		// calculate Ez field
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);

				// with pml
				E_z[cind] = gaz[cind] * D_z[cind];
				
			}
		}
		// set Ez at edges to zero for pml
		for (auto j=0; j<_mesh->reg_num_nodes_y(); j++){
			lind = _mesh->reg_inds_to_glob_ind(0, j);
			rind = _mesh->reg_inds_to_glob_ind(_mesh->reg_num_nodes_x()-1, j);
			E_z[lind] = 0.0;
			E_z[rind] = 0.0;
		}
		for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
			uind = _mesh->reg_inds_to_glob_ind(i, _mesh->reg_num_nodes_y()-1);
			dind = _mesh->reg_inds_to_glob_ind(i, 0);
			E_z[uind] = 0.0;
			E_z[dind] = 0.0;
		}
		

		// Calculate Hx
		for (auto j=0; j<_mesh->reg_num_nodes_y()-1; j++){
			for (auto i=0; i<_mesh->reg_num_nodes_x(); i++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);
				uind = _mesh->reg_inds_to_glob_ind(i, j+1);

				curle = E_z[cind] - E_z[uind];
				I_Hx[cind] = I_Hx[cind] + fi1[i]*curle;
				H_x[cind] = fj3[j]*H_x[cind]
						  + fj2[j]*_CourantFactor*(curle + I_Hx[cind]);

			}
		}
		// Calculate Hy
		for (auto j=0; j< _mesh->reg_num_nodes_y(); j++){
			for (auto i=0; i<_mesh->reg_num_nodes_x()-1; i++){
				cind = _mesh->reg_inds_to_glob_ind(i, j);
				rind = _mesh->reg_inds_to_glob_ind(i+1, j);

				curle = E_z[rind] - E_z[cind];
				I_Hy[cind] = I_Hy[cind] + fj1[j]*curle;
				H_y[cind] = fi3[i]*H_y[cind]
						  + fi2[i]*_CourantFactor*(curle + I_Hy[cind]);

			}
		}
		

		// fill in the simdata for this time step
		_simdata.add_data_at_index(n, "E_z", E_z[0]);
		//_simdata.add_data_at_index(n, "H_y", H_y[0]);

		// update iteration count
		_current_iter++;

	}


}


void FDTDSimulation::preRunCheck(){

}

void FDTDSimulation::allocate_fields(){
	cout << "allocating fields" << endl;

	D_z = new double[_mesh->nodecount()];
	I_Hx = new double[_mesh->nodecount()];
	I_Hy = new double[_mesh->nodecount()];
	E_z = new double[_mesh->nodecount()];
	H_y = new double[_mesh->nodecount()];
	H_x = new double[_mesh->nodecount()];
	gaz = new double[_mesh->nodecount()];
	//gbz = new double[_mesh.nodecount()];
	for (auto i=0; i<_mesh->nodecount(); i++){
		D_z[i] = 0.0;
		I_Hx[i] = 0.0;
		I_Hy[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
		gaz[i] = 1.0/_rel_permittivity[i];		// gaz = 1/(epsilon+(sigma*dt/eps0));
		//gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

}

void FDTDSimulation::allocate_PML(){
	cout << "allocating pml" << endl;

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
	cout << "allocating simdata" << endl;

	// set up the simulation data holder
	_simdata.bind_mesh(*_mesh);
	cout << "mesh is bound to simdata" << endl;
	_simdata.add_field("E_z");
	cout << "field is added to simdata" << endl;
	_simdata.set_time_span(0.0, _dt, _num_iters*_dt);
	cout << "allocating simdata" << endl;
	_simdata.print_summary();
}

