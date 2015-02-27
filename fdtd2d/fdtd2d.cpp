#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double eps0 = 8.854e-12;
	double c0 = 3.0e+8;
	double dx = 0.01;
	double dt = 0.5*dx/c0;

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = RegularMesh::create_regular_grid_b(dx, 0.0, 1.0, 0.0, 1.0);
	paramesh.print_summary();

	// view the mesh
	mesh_visualixer paravis;
	paravis.add_mesh(&paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	//paravis->run();


	// simple 1D simulation	
	double num_iters = 50;	

	// set up the simulation data holder
	SimulationData simdata;
	simdata.bind_mesh(paramesh);
	simdata.add_field("E_z");
	simdata.add_field("H_y");
	simdata.set_time_span(0.0, dt, num_iters*dt);
	simdata.print_summary();

	// calculate E using centered difference except for on the boundaries
	double * E_z = new double[paramesh.reg_num_nodes_x()*paramesh.reg_num_nodes_y()];
	double * H_y = new double[paramesh.reg_num_nodes_x()*paramesh.reg_num_nodes_y()];
	double * H_x = new double[paramesh.reg_num_nodes_x()*paramesh.reg_num_nodes_y()];
	for (auto i=0; i<paramesh.reg_num_nodes_x()*paramesh.reg_num_nodes_y(); i++){
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
	}

	// define some constants (for vacuum)
	double Chxh = 1.0;
	double Chxe = dt/dx;
	double Chyh = 1.0;
	double Chye = dt/dx;
	double Ceze = 1.0;
	double Cezh = dt/dx;
	
	double tcur, pulse, t0=10.0, spread = 12.0;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=1; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;

		// update E field
		for (auto i=0; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<paramesh.reg_num_nodes_y()-1; j++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);
				lind = paramesh.reg_inds_to_glob_ind(i, j-1);
				rind = paramesh.reg_inds_to_glob_ind(i, j+1);
				uind = paramesh.reg_inds_to_glob_ind(i+1, j);
				dind = paramesh.reg_inds_to_glob_ind(i-1, j);

				//H_y[cind] = H_y[cind] + 0.5*(E_x[cind] - E_x[rind]);

				E_z[cind] = Ceze * E_z[cind] + Cezh * ((H_y[cind] - H_y[dind]) - (H_x[cind] - H_x[lind]));
			}
		}

		// impose a gaussian source (hard-coded)
		pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		E_z[paramesh.reg_num_nodes_x()*paramesh.reg_num_nodes_y()/2] = pulse;

		
		// update H field
		//cout << "calculating new Electric field" << endl;
		for (auto i=1; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<paramesh.reg_num_nodes_y()-1; j++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);
				lind = paramesh.reg_inds_to_glob_ind(i, j-1);
				rind = paramesh.reg_inds_to_glob_ind(i, j+1);
				uind = paramesh.reg_inds_to_glob_ind(i+1, j);
				dind = paramesh.reg_inds_to_glob_ind(i-1, j);

				//E_x[cind] += 0.5*(H_y[lind] - H_y[cind]);
				H_x[cind] = Chxh * H_x[cind] - Chxe * (E_z[rind] - E_z[cind]);
				H_y[cind] = Chyh * H_y[cind] + Chye * (E_z[uind] - E_z[cind]);
			}
		}

		

		// fill in the simdata for this time step
		simdata.add_data_at_index(n, "E_z", E_z[0]);
		simdata.add_data_at_index(n, "H_y", H_y[0]);

		tcur += dt;
	}
	//*/

	//simdata.write_HDF5("simulation_data.h5");

	// plot the evolution of E_x
	
	for (auto i=0; i<simdata.num_time_steps(); i++){
		paravis.set_colorby(&(simdata.get_data_at_index(i, "E_z")));
		paravis.run();
	}
	
	


	delete[] E_z;
	delete[] H_y;
	return 0;
}