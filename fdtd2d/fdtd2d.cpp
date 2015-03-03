#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double eps0 = 8.854e-12;
	double c0 = 3.0e+8;
	double dx = 0.005;
	double dt = 0.5*dx/c0;

	parametric_model_2d paramodel;
	paramodel.set_model_name("DielecSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	circle c1 = circle(0.1, {0.5, 0.5}, paramodel.get_material("Dielectric"));
	circle c2 = circle(0.06, {0.5, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, paramodel.get_material("Vacuum"));
	paramesh.print_summary();

	// convert the model into a mesh
	//RegularMesh paramesh;
	//paramesh = RegularMesh::create_regular_grid_b(dx, 0.0, 1.0, 0.0, 1.0);
	//paramesh.print_summary();

	// view the mesh
	mesh_visualixer paravis;
	paravis.add_mesh(&paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	paravis.set_colorby(&paramesh.data("eps_rel"));
	paravis.run();


	// simple 1D simulation	
	double num_iters = 800;	

	// set up the simulation data holder
	SimulationData simdata;
	simdata.bind_mesh(paramesh);
	simdata.add_field("E_z");
	simdata.add_field("H_y");
	simdata.set_time_span(0.0, dt, num_iters*dt);
	simdata.print_summary();

	// calculate E using centered difference except for on the boundaries
	double * D_z = new double[paramesh.nodecount()];
	double * I_z = new double[paramesh.nodecount()];
	double * E_z = new double[paramesh.nodecount()];
	double * H_y = new double[paramesh.nodecount()];
	double * H_x = new double[paramesh.nodecount()];
	double * gaz = new double[paramesh.nodecount()];
	double * gbz = new double[paramesh.nodecount()];
	const double * epsilon_rel = &paramesh.data("eps_rel");
	for (auto i=0; i<paramesh.nodecount(); i++){
		D_z[i] = 0.0;
		I_z[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
		gaz[i] = 1.0/epsilon_rel[i];		// gaz = 1/(epsilon+(sigma*dt/eps0));
		gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

	
	double tcur, pulse, t0=10.0, spread = 12.0, srcfreq=1.0e+9;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=0; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;

		// update E field
		for (auto i=1; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<paramesh.reg_num_nodes_y()-1; j++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);
				lind = paramesh.reg_inds_to_glob_ind(i, j-1);
				rind = paramesh.reg_inds_to_glob_ind(i, j+1);
				uind = paramesh.reg_inds_to_glob_ind(i+1, j);
				dind = paramesh.reg_inds_to_glob_ind(i-1, j);

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
		E_z[paramesh.nearest_node(0.25, 0.9)] = pulse;
		E_z[paramesh.nearest_node(0.25, 0.8)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*30/180); ;
		E_z[paramesh.nearest_node(0.25, 0.7)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*60/180); ;
		E_z[paramesh.nearest_node(0.25, 0.6)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*90/180); ;
		E_z[paramesh.nearest_node(0.25, 0.5)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*120/180); ;
		E_z[paramesh.nearest_node(0.25, 0.4)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*150/180); ;
		E_z[paramesh.nearest_node(0.25, 0.3)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*180/180); ;
		E_z[paramesh.nearest_node(0.25, 0.2)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*210/180); ;
		E_z[paramesh.nearest_node(0.25, 0.1)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*240/180); ;

		
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
				H_x[cind] = H_x[cind] - 0.5 * (E_z[uind] - E_z[cind]);
				H_y[cind] = H_y[cind] + 0.5 * (E_z[rind] - E_z[cind]);
			}
		}

		// fill in the simdata for this time step
		simdata.add_data_at_index(n, "E_z", E_z[0]);
		simdata.add_data_at_index(n, "H_y", H_y[0]);

		tcur += dt;
	}
	//*/

	//simdata.write_HDF5("free_space_gaussian.h5");

	
	// visualize the simulation
	simulation_visualixer simvis;
	simvis.bind_simulation(simdata);
	simvis.set_colorby_field("E_z");
	simvis.set_color_interpolation(false);
	simvis.set_frequency_Hz(30);
	simvis.run();
	//*/
	
	

	delete[] D_z;
	delete[] I_z;
	delete[] E_z;
	delete[] H_y;
	delete[] H_x;
	delete[] gaz;
	delete[] gbz;

	return 0;
}
