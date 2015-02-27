#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double eps0 = 8.854e-12;
	double c0 = 3.0e+8;
	double dx = 0.01;
	double dt = 0.5*dx/c0;

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = RegularMesh::create_regular_grid_b(dx, 0.0, 1.0);
	paramesh.print_summary();

	// view the mesh
	mesh_visualixer paravis;
	paravis.add_mesh(&paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	//paravis->run();


	// simple 1D simulation	
	double num_iters = 300;	

	// set up the simulation data holder
	SimulationData simdata;
	simdata.bind_mesh(paramesh);
	simdata.add_field("E_x");
	simdata.add_field("H_y");
	simdata.set_time_span(0.0, dt, num_iters*dt);
	simdata.print_summary();

	// calculate E using centered difference except for on the boundaries
	double * E_x = new double[paramesh.reg_num_nodes_x()];
	double * H_y = new double[paramesh.reg_num_nodes_x()];
	H_y[0] = 0.0;
	E_x[0] = 0.0;
	H_y[paramesh.reg_num_nodes_x()-1] = 0.0;
	E_x[paramesh.reg_num_nodes_x()-1] = 0.0;


	
	double tcur, pulse, t0=40.0, spread = 12.0;
	unsigned int cind, lind, rind;
	for (auto n=1; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;

		// impose boundary conditions
		//H_y[paramesh.reg_num_nodes_x()-1] = H_y[paramesh.reg_num_nodes_x()-2];
		H_y[0] = H_y[1];

		// update H field
		for (auto i=0; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
				cind = paramesh.reg_inds_to_glob_ind(i);
				rind = paramesh.reg_inds_to_glob_ind(i+1);

				H_y[cind] = H_y[cind] + 0.5*(E_x[cind] - E_x[rind]);

				//E_x[cind] += 0.5*(H_y[lind] - H_y[cind]);
		}

		// impose boundary conditions
		E_x[0] = E_x[1];
		//E_x[paramesh.reg_num_nodes_x()-1] = E_x[paramesh.reg_num_nodes_x()-2];

		// update E field
		//cout << "calculating new Electric field" << endl;
		for (auto i=1; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
				cind = paramesh.reg_inds_to_glob_ind(i);
				lind = paramesh.reg_inds_to_glob_ind(i-1);

				E_x[cind] += 0.5*(H_y[lind] - H_y[cind]);
		}

		// impose a gaussian source
		pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		E_x[paramesh.reg_num_nodes_x()/2] += pulse;

		

		// fill in the simdata for this time step
		simdata.add_data_at_index(n, "E_x", E_x[0]);
		simdata.add_data_at_index(n, "H_y", H_y[0]);

		tcur += dt;
	}
	//*/

	simdata.write_HDF5("outfile_test.h5");

	// plot the evolution of E_x
	simulation_visualixer simvis;
	simvis.bind_simulation(simdata);
	simvis.set_colorby_field("E_x");
	simvis.set_frequency_Hz(100);
	simvis.run();
	//cout << "finished running..." << endl;
	/*
	for (auto i=0; i<simdata.num_time_steps(); i++){
		cout << "on time step " << i << "/" << simdata.num_time_steps() << "\r" << flush;
		paravis.set_colorby(&(simdata.get_data_at_index(i, "E_x")));
		paravis.run();
	}
	*/
	


	delete[] E_x;
	delete[] H_y;
	return 0;
}