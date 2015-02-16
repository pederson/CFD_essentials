//#include "../Visualixer.hpp" 
#include "../VisualixerMesh.hpp"
#include "../GeometricObject.hpp"
#include "../Mesh.hpp" 
#include "../Converter.hpp" 
//#include "EquationTerm.hpp"
//#include "Equation.hpp"
//#include "Simulation.hpp"
#include <petscksp.h>
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
	
	
	// loop over the internal nodes adding the laplace operator to the matrix
	PetscScalar * mat = new PetscScalar [paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) mat[i] = 0;

	unsigned int cind, lind, rind, uind, dind;
	for (auto i=1; i<paramesh->reg_num_nodes_x()-1; i++){ // columns
		for (auto j=1; j<paramesh->reg_num_nodes_y()-1; j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			lind = paramesh->reg_inds_to_glob_ind(i-1,j);
			rind = paramesh->reg_inds_to_glob_ind(i+1,j);
			uind = paramesh->reg_inds_to_glob_ind(i,j+1);
			dind = paramesh->reg_inds_to_glob_ind(i,j-1);
			mat[cind*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() + cind] = -4.0;
			mat[cind*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() + lind] = 1.0;
			mat[cind*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() + rind] = 1.0;
			mat[cind*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() + uind] = 1.0;
			mat[cind*paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() + dind] = 1.0;
		}
	}



	// construct the right side using the density field
	double q_electron = -1.6e-19, eps0 = 8.854e-12, m_electron = 9.11e-31;
	paramesh->print_summary();
	const double * reldata = &paramesh->data("e_density");
	double * rhs = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_x(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_y(); j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			//cout << "I: " << i << " J: " << j << " density: " << reldata[cind] << endl;
			rhs[cind] = -reldata[cind]*q_electron/eps0*dx*dx;
		}
	}

	
	

	// now solve the system (PETSc KSP)
	Vec x,b;   /* approx solution, rhs */
    Mat A;     /* linear system matrix */
    KSP ksp;   /* linear solver context */
    PC pc;     /* preconditioner context */
    int size;

	cout << "starting petsc stuff" << endl;
	MPI_Init(&argc,&argv); // this is required before PetscInitialize
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    cout << "initialized petsc" << endl;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    cout << "initialized mpi comm with " << size << " processes" << endl;

    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x());
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
    

    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(),paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x());
    MatSetFromOptions(A);
    MatSetUp(A);

    
    cout << "about to init matrix" << endl;

	// set the matrix values for the laplace operator
	PetscScalar oper[5];
	oper[0] = -1.0; oper[1] = -1.0; oper[2] = 4.0; oper[3] = -1.0; oper[4] = -1.0;
	PetscInt cols[5], row;
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
			MatSetValues(A, 1, &row, 5, cols, oper, INSERT_VALUES);
		}
	}
	cout << endl;


	/*
	// insert boundary conditions
	PetscInt cl;
	PetscScalar bnds[4];
	bnds[0] = 0.0; bnds[1] = 0.0; bnds[2] = 1.0; bnds[3] = 0.0;
	PetscInt cls[4], rw;
	for (auto j=1; j<paramesh->reg_num_nodes_y()-1; j++){ // rows
		cind = paramesh->reg_inds_to_glob_ind(paramesh->reg_num_nodes_x()-1, j);
		lind = paramesh->reg_inds_to_glob_ind(paramesh->reg_num_nodes_x()-2, j);
		uind = paramesh->reg_inds_to_glob_ind(paramesh->reg_num_nodes_x()-1, j+1);
		dind = paramesh->reg_inds_to_glob_ind(paramesh->reg_num_nodes_x()-1, j-1);
		cls[0] = lind;
		cls[1] = dind;
		cls[2] = cind;
		cls[3] = uind;

		//rw = j*paramesh->reg_num_nodes_y() + paramesh->reg_num_nodes_x();
		rw = cind;
		MatSetValues(A, 1, &rw, 4, cls, bnds, INSERT_VALUES);

	}
	// write in the boundary conditions
	// -100 V on the right side of the domain
	for (auto j=0; j<paramesh->reg_num_nodes_y(); j++){ // rows
		cind = paramesh->reg_inds_to_glob_ind(paramesh->reg_num_nodes_x()-1,j);
		//cout << "I: " << i << " J: " << j << endl;
		rhs[cind] = -100.0;

	}
	//*/

	for(auto i=0;i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x();i++) VecSetValues(b,1,&i,&rhs[i],INSERT_VALUES);
    /* need to assemble after setting values! do necessary
       message passing etc to propagate matrix to all ranks */
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    cout << "assembled rhs" << endl;
	
	
	
    /* need to assemble matrix for the same reasons as above */
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    cout << "assembled matrix" << endl;

	KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCJACOBI);
    KSPSetTolerances(ksp,1e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp,b,x);
    cout << "solved equation" << endl;

    //MatView(A,PETSC_VIEWER_STDOUT_WORLD);
    MatView(A,PETSC_VIEWER_DRAW_WORLD); // this will draw the non-zero entries of the matrix
    cout << "enter something: " << endl;
    string input = "";
    //getline(cin, input);
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);
    //VecView(x,PETSC_VIEWER_STDOUT_WORLD);
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

    PetscInt sz;
    VecGetSize(x, &sz);
    cout << "NUM in vector solution: " << sz << "    NUM in rhs: " << paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x() << endl;

    // collect the value from the solver
    double * soln = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
    PetscInt * inds = new PetscInt[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
    for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) inds[i] = i;
    cout << "gonna do some ish" << endl;
    VecGetValues(x, paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(), inds, soln);
    cout << "did some ish" << endl;

    for (auto i=sz-10; i<sz; i++) cout << "Solution[" << i << "]: " << soln[i] << endl;
    cout << "finished that" << endl;

	// visualize the solution by putting it in the mesh
	paramesh->add_phys_property("potential", soln);
	cout << "added property" << endl;
	paramesh->print_summary();
	paravis->set_colorby(&paramesh->data("potential"));
	cout << "set the colorby" << endl;
	paravis->run();
	cout << "finished running" << endl;



	/*
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

			E_x[cind] = (soln[rind] - soln[lind])/(2*dx);
			E_y[cind] = (soln[uind] - soln[dind])/(2*dx);
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

	double dt = 1e-9; // one nanosecond time step
	double num_iters = 5;
	// use forward euler to advance velocity in time
	// NEED TO FIX THIS TO DEAL WITH BOUNDARIES WHEN CENTRAL DIFFERENCE OPERATOR IS USED
	
	for (auto n=0; n<num_iters; n++){

		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				lind = paramesh->reg_inds_to_glob_ind(j, i-1);
				rind = paramesh->reg_inds_to_glob_ind(j, i+1);
				uind = paramesh->reg_inds_to_glob_ind(j+1, i);
				dind = paramesh->reg_inds_to_glob_ind(j-1, i);

				// calculate updated velocity
				vel_x[cind] = vel_x_old[cind] + dt*(-q_electron*E_x[cind]/m_electron - 760*5.3e+9*vel_x_old[cind]);
				vel_x_old[cind] = vel_x[cind];

				vel_y[cind] = vel_y_old[cind] + dt*(-q_electron*E_y[cind]/m_electron - 760*5.3e+9*vel_y_old[cind]);
				vel_y_old[cind] = vel_y[cind];
			}
		}

		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				lind = paramesh->reg_inds_to_glob_ind(j, i-1);
				rind = paramesh->reg_inds_to_glob_ind(j, i+1);
				uind = paramesh->reg_inds_to_glob_ind(j+1, i);
				dind = paramesh->reg_inds_to_glob_ind(j-1, i);

				// calculate updated density field
				n_e[cind] = n_e_old[cind] + dt * (n_e_old[cind]*(vel_x[rind]-vel_x[lind])/(2*dx)
										  + vel_x[cind]*(n_e_old[rind]-n_e_old[lind])/(2*dx)
										  + n_e_old[cind]*(vel_y[uind]-vel_y[dind])/(2*dx)
										  + vel_y[cind]*(n_e_old[uind]-n_e_old[dind])/(2*dx)
										  );
				n_e_old[cind] = n_e[cind]; // IS THIS A PROBLEM IF I UPDATE HERE?
			}
		}

		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
			for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
				cind = paramesh->reg_inds_to_glob_ind(j, i);
				rhs[cind] = -n_e[cind]*q_electron/eps0*dx*dx;
			}
		}
		// solve poisson's equation
		for(auto i=0;i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x();i++) VecSetValues(b,1,&i,&rhs[i],INSERT_VALUES);
	    // need to assemble after setting values! do necessary
	    //   message passing etc to propagate matrix to all ranks 
	    VecAssemblyBegin(b);
	    VecAssemblyEnd(b);
	    KSPSolve(ksp,b,x);


		// update E field
		for (auto j=1; j<paramesh->reg_num_nodes_x()-1; j++){ // cols
		for (auto i=1; i<paramesh->reg_num_nodes_y()-1; i++){
			cind = paramesh->reg_inds_to_glob_ind(j, i);
			lind = paramesh->reg_inds_to_glob_ind(j, i-1);
			rind = paramesh->reg_inds_to_glob_ind(j, i+1);
			uind = paramesh->reg_inds_to_glob_ind(j+1, i);
			dind = paramesh->reg_inds_to_glob_ind(j-1, i);

			E_x[cind] = (soln[rind] - soln[lind])/(2*dx);
			E_y[cind] = (soln[uind] - soln[dind])/(2*dx);
		}
	}
	}
	*/



	/* this is what I would do if I had all the classes ready
	// create an equation
	// del^2(potential) = -e*(e_density)/epsilon
	//Equation poiss = Equation(EQUATION_POISSON);

	// set boundary conditions


	// create a simulation

	// run the simulation

	// visualize the results
	*/



	// cleanup
    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    PetscFinalize();

	delete[] rhs;
	delete[] mat;
	delete paravis;


	return 0;
}