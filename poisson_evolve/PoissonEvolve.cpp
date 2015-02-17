//#include "../Visualixer.hpp" 
#include "../VisualixerMesh.hpp"
#include "../GeometricObject.hpp"
#include "../Mesh.hpp" 
#include "../Converter.hpp" 
//#include "EquationTerm.hpp"
//#include "Equation.hpp"
//#include "Simulation.hpp"
#include "../LinAlgWrapper.hpp"
//#include <petscksp.h>
#include <omp.h>

using namespace std;

int main(int argc, char * argv[]){

	// define the geometry model that will be used
	parametric_model_2d my_param2;
	my_param2.set_model_name("GaussianDistElectrons");
	my_param2.add_physical_property("e_density");
	my_param2.add_material("Air", {1.0e+10});
	my_param2.add_material("Dielectric", {5.0e+13});
	my_param2.add_material("Dielectric2", {9.0});
	//gaussian_2d ga1 = gaussian_2d(0.3, 0.3, 1000.0, 1.0e+16, {0.0, 0.0}); // electron density in #/m^3
	circle c1 = circle(0.03, {0.0, 0.04}, my_param2.get_material("Dielectric"));
	circle c2 = circle(0.03, {0.0, -0.04}, my_param2.get_material("Dielectric"));
	//my_param2.add_object(&ga1);
	my_param2.add_object(&c1);
	my_param2.add_object(&c2);

	// convert the model into a mesh
	Static_Mesh * paramesh;
	double dx = 0.002;
	paramesh = build_simple_mesh_2d(&my_param2, dx, -0.1, 0.1, -0.15, 0.15, my_param2.get_material("Air"));
	
	// view the mesh
	mesh_visualixer * paravis = new mesh_visualixer();
	paravis->add_mesh(paramesh);
	paravis->set_color_ramp(CRamp::DIVERGENT_9);
	paravis->set_colorby(&paramesh->data("e_density"));
	paravis->run();
	//delete paramesh; // WHY DOES THIS COMPILE IF THIS IS UNCOMMENTED???


	// this is a testing section
		// solve the poisson equation using a finite difference method
	//Simulation mysim = Simulation(FINITE_DIFFERENCE, paramesh);

	// construct the right side using the density field
	double q_electron = -1.6e-19, eps0 = 8.854e-12, m_electron = 9.11e-31;
	paramesh->print_summary();
	unsigned int cind;
	const double * reldata = &paramesh->data("e_density");
	double * rhs = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_x(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_y(); j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			//cout << "I: " << i << " J: " << j << " density: " << reldata[cind] << endl;
			rhs[cind] = -reldata[cind]*q_electron/eps0*dx*dx;
		}
	}

	// now solve the system (linear algebra wrapper)
	LinVector rhsv(paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x());   // approx solution, rhs 
    LinMatrix A(paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(), paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x());     // linear system matrix 
    LinSolver ksp;   // linear solver context 
    int size;

    

    cout << "about to init matrix" << endl;

	// set the matrix values for the laplace operator
	unsigned int lind, rind, uind, dind;
	double oper[5];
	oper[0] = -1.0; oper[1] = -1.0; oper[2] = 4.0; oper[3] = -1.0; oper[4] = -1.0;
	int cols[5], row;
	// for each row
	for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
		for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(j, i);
			lind = paramesh->reg_inds_to_glob_ind(j, i-1);
			rind = paramesh->reg_inds_to_glob_ind(j, i+1);
			uind = paramesh->reg_inds_to_glob_ind(j+1, i);
			dind = paramesh->reg_inds_to_glob_ind(j-1, i);
			cols[0] = lind;
			cols[1] = dind;
			cols[2] = cind;
			cols[3] = uind;
			cols[4] = rind;

			if (j%10 == 0) cout << "on row: " << j << " / " << paramesh->reg_num_nodes_x() << " \r" << flush;
			//cout << "yer" << endl;
			//row = j*paramesh->reg_num_nodes_x() + i; // row index in the overall matrix
			row = cind;
			A.insert_values(1, &row, 5, cols, oper);
		}
	}
	cout << endl;


	for(auto i=0;i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x();i++){
		rhsv.insert_values(1,&rhs[i], &i);
	}
    cout << "assembled rhs" << endl;
	
	LinVector x = ksp.solve(A, rhsv);
    cout << "solved equation" << endl;
    const double * dat;
    dat = &x.data();
    for (auto i=0; i<100; i++){
    	cout << "rhs[" << i << "]: " << rhsv[i] << "    " ;
    	cout << "solution[" << i << "]: " << x[i] << "    ";
    	cout << "solution[" << i << "]: " << dat[i] << endl;
    } 

    //A.draw(); // a looks fine... I checked already

    
    unsigned int sz = x.length();
    cout << "NUM in vector solution: " << sz << "    NUM in rhs: " << paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() << endl;

    // collect the value from the solver

	// visualize the solution by putting it in the mesh
	//paramesh->add_phys_property("potential", &x.data());
	cout << "added property" << endl;
	//paramesh->print_summary();
	//paravis->set_colorby(&paramesh->data("potential"));
	paravis->set_colorby(&x.data());
	cout << "set the colorby" << endl;
	paravis->run();
	cout << "finished running" << endl;
	
	

	
	//************* UNSTEADY PROBLEM ****************
	// using dv_e/dt = -eE/m_e - nu_c*v_e
	// and   nu_c = 5.3e+9 * P_atm (in torr)
	// 
	// P_atm = 760 torr
	//
	// where E = -grad(phi)
	//
	// then solve the continuity equation
	// dn_e/dt + del . (n_e*v_e) = 0

	// calculate E using centered difference except for on the boundaries
	double * E_x = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	double * E_y = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
		for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
			cind = paramesh->reg_inds_to_glob_ind(j, i);
			lind = paramesh->reg_inds_to_glob_ind(j, i-1);
			rind = paramesh->reg_inds_to_glob_ind(j, i+1);
			uind = paramesh->reg_inds_to_glob_ind(j+1, i);
			dind = paramesh->reg_inds_to_glob_ind(j-1, i);

			E_x[cind] = (x[rind] - x[lind])/(2*dx);
			E_y[cind] = (x[uind] - x[dind])/(2*dx);
		}
	}
	// deal with the boundaries by doing a forward (?) difference

	double * vel_x = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	double * vel_y = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	double * vel_x_old = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	double * vel_y_old = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	// initially, the velocities are 0 (the plasma is stationary)
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++){
		vel_x[i] = 0.0;
		vel_y[i] = 0.0;
		vel_x_old[i] = 0.0;
		vel_y_old[i] = 0.0;
	}

	double * n_e = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	double * n_e_old = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) n_e_old[i] = reldata[i];

	double dt = 1e-11; // one nanosecond time step
	double num_iters = 5;
	// use forward euler to advance velocity in time
	// NEED TO FIX THIS TO DEAL WITH BOUNDARIES WHEN CENTRAL DIFFERENCE OPERATOR IS USED
	
	for (auto n=0; n<num_iters; n++){

		cout << "calculating new velocity" << endl;
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				lind = paramesh->reg_inds_to_glob_ind(j, i-1);
				rind = paramesh->reg_inds_to_glob_ind(j, i+1);
				uind = paramesh->reg_inds_to_glob_ind(j+1, i);
				dind = paramesh->reg_inds_to_glob_ind(j-1, i);

				// calculate updated velocity
				vel_x[cind] = vel_x_old[cind] + dt*(-q_electron*E_x[cind]/m_electron - 760*5.3e+9*vel_x_old[cind]);
				//vel_x_old[cind] = vel_x[cind];

				vel_y[cind] = vel_y_old[cind] + dt*(-q_electron*E_y[cind]/m_electron - 760*5.3e+9*vel_y_old[cind]);
				//vel_y_old[cind] = vel_y[cind];
			}
		}
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				vel_x_old[cind] = vel_x[cind];
				vel_y_old[cind] = vel_y[cind];
			}
		}

		cout << "plotting vel_x" << endl;
		paravis->set_colorby(vel_x);
		paravis->run();
		cout << "plotting vel_y" << endl;
		paravis->set_colorby(vel_y);
		paravis->run();

		cout << "calculating new density" << endl;
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				lind = paramesh->reg_inds_to_glob_ind(j, i-1);
				rind = paramesh->reg_inds_to_glob_ind(j, i+1);
				uind = paramesh->reg_inds_to_glob_ind(j+1, i);
				dind = paramesh->reg_inds_to_glob_ind(j-1, i);

				// calculate updated density field
				n_e[cind] = n_e_old[cind] + dt * 
										  (
										  	n_e_old[cind]*(vel_x[rind]-vel_x[lind])/(2*dx)
										  + vel_x[cind]*(n_e_old[rind]-n_e_old[lind])/(2*dx)
										  + n_e_old[cind]*(vel_y[uind]-vel_y[dind])/(2*dx)
										  + vel_y[cind]*(n_e_old[uind]-n_e_old[dind])/(2*dx)
										  );
				//n_e_old[cind] = n_e[cind]; // IS THIS A PROBLEM IF I UPDATE HERE?
			}
		}
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				n_e_old[cind] = n_e[cind];
			}
		}

		cout << "plotting electron density" << endl;
		paravis->set_colorby(n_e);
		paravis->run();

		cout << "calculating new rhs and solving for potential" << endl;
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				rhs[cind] = -n_e[cind]*q_electron/eps0*dx*dx;
			}
		}
		// solve poisson's equation
		for(auto i=0;i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x();i++) rhsv[i] = rhs[i];
	    // need to assemble after setting values! do necessary
	    //   message passing etc to propagate matrix to all ranks 
	    x = ksp.solve(A, rhsv);


		// update E field
		cout << "calculating new Electric field" << endl;
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				lind = paramesh->reg_inds_to_glob_ind(j, i-1);
				rind = paramesh->reg_inds_to_glob_ind(j, i+1);
				uind = paramesh->reg_inds_to_glob_ind(j+1, i);
				dind = paramesh->reg_inds_to_glob_ind(j-1, i);

				E_x[cind] = (x[rind] - x[lind])/(2*dx);
				E_y[cind] = (x[uind] - x[dind])/(2*dx);
			}
		}
		cout << "plotting E_x" << endl;
		paravis->set_colorby(E_x);
		paravis->run();
		cout << "plotting E_y" << endl;
		paravis->set_colorby(E_y);
		paravis->run();


	}
	//*/



	/* this is what I would do if I had all the classes ready
	// create an equation
	// del^2(potential) = -e*(e_density)/epsilon
	//Equation poiss = Equation(EQUATION_POISSON);

	// set boundary conditions


	// create a simulation

	// run the simulation

	// visualize the results
	*/

	clearLinAlg();
	delete[] rhs;
	delete[] vel_x;
	delete[] vel_y;
	delete[] vel_x_old;
	delete[] vel_y_old;
	delete[] n_e;
	delete[] n_e_old;
	delete[] E_x;
	delete[] E_y;
	delete paravis;
	return 0;
}