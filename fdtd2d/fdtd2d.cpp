#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	const double eps0 = 8.854e-12;
	double c0 = 3.0e+8;
	double dx = 0.005;
	double dt = 0.5*dx/c0;

	// empty space model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("EmptySpace");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	circle c1 = circle(0.1, {0.5, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&c1);
	//*/

	// ring model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("DielecSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	circle c1 = circle(0.1, {0.5, 0.5}, paramodel.get_material("Dielectric"));
	circle c2 = circle(0.06, {0.5, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	//*/
	

	// holey waveguide model
	
	parametric_model_2d paramodel;
	paramodel.set_model_name("DielecSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	rectangle r1 = rectangle(1.0, 0.3, {0.5, 0.5}, paramodel.get_material("Dielectric"));
	circle c1 = circle(0.06, {0.2, 0.5}, paramodel.get_material("Vacuum"));
	circle c2 = circle(0.06, {0.4, 0.5}, paramodel.get_material("Vacuum"));
	circle c3 = circle(0.06, {0.6, 0.5}, paramodel.get_material("Vacuum"));
	circle c4 = circle(0.06, {0.8, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	//*/


	// holey waveguide defect model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("DielecSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	rectangle r1 = rectangle(1.0, 0.3, {0.5, 0.5}, paramodel.get_material("Dielectric"));
	circle c1 = circle(0.05, {0.15, 0.5}, paramodel.get_material("Vacuum"));
	circle c2 = circle(0.05, {0.3, 0.5}, paramodel.get_material("Vacuum"));
	circle c3 = circle(0.05, {0.7, 0.5}, paramodel.get_material("Vacuum"));
	circle c4 = circle(0.05, {0.85, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	//*/
	

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, paramodel.get_material("Vacuum"));
	paramesh.print_summary();


	// view the mesh
	mesh_visualixer paravis;
	paravis.bind_mesh(paramesh);
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
	double * I_Hx = new double[paramesh.nodecount()];
	double * I_Hy = new double[paramesh.nodecount()];
	double * E_z = new double[paramesh.nodecount()];
	double * H_y = new double[paramesh.nodecount()];
	double * H_x = new double[paramesh.nodecount()];
	double * gaz = new double[paramesh.nodecount()];
	double * gbz = new double[paramesh.nodecount()];
	const double * epsilon_rel = &paramesh.data("eps_rel");
	for (auto i=0; i<paramesh.nodecount(); i++){
		D_z[i] = 0.0;
		I_z[i] = 0.0;
		I_Hx[i] = 0.0;
		I_Hy[i] = 0.0;
		H_x[i] = 0.0;
		H_y[i] = 0.0;
		E_z[i] = 0.0;
		gaz[i] = 1.0/epsilon_rel[i];		// gaz = 1/(epsilon+(sigma*dt/eps0));
		gbz[i] = 0.0;		// gbz = sigma*dt/eps0;
	}

	// PML related
	float npml, xnum, xd, xxn, xn;
	npml = 8;
	double * gi2 = new double[paramesh.reg_num_nodes_x()];
	double * gi3 = new double[paramesh.reg_num_nodes_x()];
	double * gj2 = new double[paramesh.reg_num_nodes_y()];
	double * gj3 = new double[paramesh.reg_num_nodes_y()];
	double * fi1 = new double[paramesh.reg_num_nodes_x()];
	double * fi2 = new double[paramesh.reg_num_nodes_x()];
	double * fi3 = new double[paramesh.reg_num_nodes_x()];
	double * fj1 = new double[paramesh.reg_num_nodes_y()];
	double * fj2 = new double[paramesh.reg_num_nodes_y()];
	double * fj3 = new double[paramesh.reg_num_nodes_y()];
	for (auto i=0; i<paramesh.reg_num_nodes_x(); i++){
		gi2[i] = 1.0;
		gi3[i] = 1.0;
		fi1[i] = 0.0;
		fi2[i] = 1.0;
		fi3[i] = 1.0;
	}
	for (auto i=0; i<paramesh.reg_num_nodes_y(); i++){
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
		gi2[paramesh.reg_num_nodes_x()-1-i] = 1.0/(1.0+xn);
		gi3[paramesh.reg_num_nodes_x()-1-i] = (1.0-xn)/(1.0+xn);

		xxn = (xnum-0.5)/xd;
		xn = 0.25*xxn*xxn*xxn;
		fi1[i] = xn;
		fi2[i] = 1.0/(1.0+xn);
		fi3[i] = (1.0-xn)/(1.0+xn);
		fi1[paramesh.reg_num_nodes_x()-2-i] = xn;
		fi2[paramesh.reg_num_nodes_x()-2-i] = 1.0/(1.0+xn);
		fi3[paramesh.reg_num_nodes_x()-2-i] = (1.0-xn)/(1.0+xn);


		// top and bottom boundaries
		xxn = xnum/xd;
		xn = 0.33*xxn*xxn*xxn;
		gj2[i] = 1.0/(1.0+xn);
		gj3[i] = (1.0-xn)/(1.0+xn);
		gj2[paramesh.reg_num_nodes_y()-1-i] = 1.0/(1.0+xn);
		gj3[paramesh.reg_num_nodes_y()-1-i] = (1.0-xn)/(1.0+xn);

		xxn = (xnum-0.5)/xd;
		xn = 0.25*xxn*xxn*xxn;
		fj1[i] = xn;
		fj2[i] = 1.0/(1.0+xn);
		fj3[i] = (1.0-xn)/(1.0+xn);
		fj1[paramesh.reg_num_nodes_y()-2-i] = xn;
		fj2[paramesh.reg_num_nodes_y()-2-i] = 1.0/(1.0+xn);
		fj3[paramesh.reg_num_nodes_y()-2-i] = (1.0-xn)/(1.0+xn);
	}

	double CourantFactor = 0.5;
	double tcur, pulse, t0=10.0, spread = 12.0, srcfreq=1.0e+9, curle;
	unsigned int cind, lind, rind, uind, dind ;
	for (auto n=0; n<num_iters; n++){
		tcur = dt;

		cout << "on time step " << n << "/" << num_iters-1 << "\r" << flush;


		// update D field
		for (auto i=1; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<paramesh.reg_num_nodes_y()-1; j++){

				cind = paramesh.reg_inds_to_glob_ind(i, j);
				lind = paramesh.reg_inds_to_glob_ind(i-1, j);
				rind = paramesh.reg_inds_to_glob_ind(i+1, j);
				uind = paramesh.reg_inds_to_glob_ind(i, j+1);
				dind = paramesh.reg_inds_to_glob_ind(i, j-1);

				/* without pml
				D_z[cind] = D_z[cind] 
						  + 0.5 * (H_y[cind] - H_y[lind] - H_x[cind] + H_x[dind]);
				E_z[cind] = gaz[cind] *(D_z[cind] - I_z[cind]);
				I_z[cind] = I_z[cind] + gbz[cind] * E_z[cind];
				*/

				
				D_z[cind] = gi3[i]*gj3[j]*D_z[cind]
						  + gi2[i]*gj2[j]*CourantFactor*(H_y[cind] - H_y[lind] - H_x[cind] + H_x[dind]);
				//E_z[cind] = gaz[cind] *(D_z[cind] - I_z[cind]);
				//I_z[cind] = I_z[cind] + gbz[cind] * E_z[cind];
				
			}
		}



		// impose a gaussian source (hard-coded)
		//pulse = exp(-0.5*(t0-n)*(t0-n)/spread/spread);
		pulse = sin(2*VX_PI*srcfreq*n*dt); // oscillatory source
		D_z[paramesh.nearest_node(0.05, 0.5)] = pulse;
		//E_z[paramesh.nearest_node(0.5, 0.5)] = pulse; without pml
		/*
		E_z[paramesh.nearest_node(0.25, 0.9)] = pulse;
		E_z[paramesh.nearest_node(0.25, 0.8)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*30/180); ;
		E_z[paramesh.nearest_node(0.25, 0.7)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*60/180); ;
		E_z[paramesh.nearest_node(0.25, 0.6)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*90/180); ;
		E_z[paramesh.nearest_node(0.25, 0.5)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*120/180); ;
		E_z[paramesh.nearest_node(0.25, 0.4)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*150/180); ;
		E_z[paramesh.nearest_node(0.25, 0.3)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*180/180); ;
		E_z[paramesh.nearest_node(0.25, 0.2)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*210/180); ;
		E_z[paramesh.nearest_node(0.25, 0.1)] = sin(2*VX_PI*srcfreq*n*dt + VX_PI*240/180); ;
		*/

		
		// calculate Ez field
		for (auto i=0; i<paramesh.reg_num_nodes_x(); i++){ // cols
			for (auto j=0; j<paramesh.reg_num_nodes_y(); j++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);

				// with pml
				E_z[cind] = gaz[cind] * D_z[cind];
				
			}
		}
		// set Ez at edges to zero for pml
		for (auto j=0; j<paramesh.reg_num_nodes_y(); j++){
			lind = paramesh.reg_inds_to_glob_ind(0, j);
			rind = paramesh.reg_inds_to_glob_ind(paramesh.reg_num_nodes_x()-1, j);
			E_z[lind] = 0.0;
			E_z[rind] = 0.0;
		}
		for (auto i=0; i<paramesh.reg_num_nodes_x(); i++){
			uind = paramesh.reg_inds_to_glob_ind(i, paramesh.reg_num_nodes_y()-1);
			dind = paramesh.reg_inds_to_glob_ind(i, 0);
			E_z[uind] = 0.0;
			E_z[dind] = 0.0;
		}
		
		
		/*
		// update H field
		//cout << "calculating new Electric field" << endl;
		for (auto i=1; i<paramesh.reg_num_nodes_x()-1; i++){ // cols
			for (auto j=1; j<paramesh.reg_num_nodes_y()-1; j++){

				cind = paramesh.reg_inds_to_glob_ind(i, j);
				lind = paramesh.reg_inds_to_glob_ind(i-1, j);
				rind = paramesh.reg_inds_to_glob_ind(i+1, j);
				uind = paramesh.reg_inds_to_glob_ind(i, j+1);
				dind = paramesh.reg_inds_to_glob_ind(i, j-1);

				H_x[cind] = H_x[cind] - 0.5 * (E_z[uind] - E_z[cind]);
				H_y[cind] = H_y[cind] + 0.5 * (E_z[rind] - E_z[cind]);
				
			}
		}
		*/

		
		// Calculate Hx
		for (auto j=0; j<paramesh.reg_num_nodes_y()-1; j++){
			for (auto i=0; i<paramesh.reg_num_nodes_x(); i++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);
				uind = paramesh.reg_inds_to_glob_ind(i, j+1);

				curle = E_z[cind] - E_z[uind];
				I_Hx[cind] = I_Hx[cind] + fi1[i]*curle;
				H_x[cind] = fj3[j]*H_x[cind]
						  + fj2[j]*CourantFactor*(curle + I_Hx[cind]);

			}
		}
		// Calculate Hy
		for (auto j=0; j< paramesh.reg_num_nodes_y(); j++){
			for (auto i=0; i<paramesh.reg_num_nodes_x()-1; i++){
				cind = paramesh.reg_inds_to_glob_ind(i, j);
				rind = paramesh.reg_inds_to_glob_ind(i+1, j);

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
	//*/

	//simdata.write_HDF5("freespace_pml.h5");

	
	// visualize the simulation
	simulation_visualixer simvis;
	simvis.bind_simulation(simdata);
	simvis.set_colorby_field("E_z");
	simvis.set_color_interpolation(false);
	simvis.set_color_alpha(&paramesh.data("eps_rel"));
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
