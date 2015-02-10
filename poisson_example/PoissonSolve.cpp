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
	my_param2.add_material("Air", {1.0});
	my_param2.add_material("Dielectric", {5.0});
	my_param2.add_material("Dielectric2", {9.0});
	gaussian_2d ga1 = gaussian_2d(0.3, 0.3, 1000.0, 1.0e+16, {0.0, 0.0}); // electron density in #/m^3
	my_param2.add_object(&ga1);

	// convert the model into a mesh
	Static_Mesh * paramesh;
	paramesh = build_simple_mesh_2d(&my_param2, 0.01, -1.0, 1.0, -1.0, 1.0, my_param2.get_material("Air"));
	
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
	double ** mat = new double *[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) mat[i] = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); j++) mat[i][j] = 0;
	}

	unsigned int cind, lind, rind, uind, dind;
	for (auto i=1; i<paramesh->reg_num_nodes_x()-1; i++){ // columns
		for (auto j=1; j<paramesh->reg_num_nodes_y()-1; j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			lind = paramesh->reg_inds_to_glob_ind(i-1,j);
			rind = paramesh->reg_inds_to_glob_ind(i+1,j);
			uind = paramesh->reg_inds_to_glob_ind(i,j+1);
			dind = paramesh->reg_inds_to_glob_ind(i,j-1);
			mat[cind][cind] = -4.0;
			mat[cind][lind] = 1.0;
			mat[cind][rind] = 1.0;
			mat[cind][uind] = 1.0;
			mat[cind][dind] = 1.0;
		}
	}

	// write in the boundary conditions
	// zeros for now...just testing the code

	// construct the right side using the density field
	double q_electron = -1.6e-19, eps0 = 8.854e-12;
	paramesh->print_summary();
	const double * reldata = &paramesh->data("e_density");
	double * rhs = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_x(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_y(); j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			//cout << "I: " << i << " J: " << j << endl;
			rhs[cind] = -reldata[cind]*q_electron/eps0;
		}
	}

	cout << "RHS[0]: " << rhs[0] << endl;
	
	
	

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
    for(auto i=0;i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x();i++) VecSetValues(b,1,&i,&rhs[i],INSERT_VALUES);
    /* need to assemble after setting values! do necessary
       message passing etc to propagate matrix to all ranks */
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    cout << "assembled rhs" << endl;

    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(),paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x());
    MatSetFromOptions(A);
    MatSetUp(A);

    

    PetscInt idxm[paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y()], idxn[paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y()];;
    for (auto i=0; i<paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(); i++){
    	idxm[i] = i;
    	idxn[i] = i;
    }
    MatSetValues(A, paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(), idxm, paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(), idxn, mat[0], INSERT_VALUES);
//#pragma omp parallel for collapse(2)
    /*
    for (auto i=0; i<paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(); j++){ // rows
			MatSetValue(A,i,j,mat[i][j],INSERT_VALUES);
			//MatSetValues(A,paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y(),rowinds,1,&j,&mat[0][j],INSERT_VALUES);
			if (i%10 == 0) cout << "on row: " << i << " / " << paramesh->reg_num_nodes_x()*paramesh->reg_num_nodes_y() << " \r" << flush;
		}
	}
	*/
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
    //MatView(A,PETSC_VIEWER_DRAW_WORLD); // this will draw the non-zero entries of the matrix
    //cout << "enter something: " << endl;
    //string input = "";
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
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) delete[] mat[i];
	delete[] mat;
	delete paravis;


	return 0;
}