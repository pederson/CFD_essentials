//#include "../Visualixer.hpp" 
#include "../VisualixerMesh.hpp"
#include "../GeometricObject.hpp"
#include "../Mesh.hpp" 
#include "../Converter.hpp" 
//#include "EquationTerm.hpp"
//#include "Equation.hpp"
//#include "Simulation.hpp"
#include <petsc.h>

using namespace std;

int main(int argc, char * argv[]){

	// define the geometry model that will be used
	parametric_model_2d my_param2;
	my_param2.set_model_name("GaussianDistElectrons");
	my_param2.add_physical_property("e_density");
	my_param2.add_material("Air", {1.0});
	my_param2.add_material("Dielectric", {5.0});
	my_param2.add_material("Dielectric2", {9.0});
	gaussian_2d ga1 = gaussian_2d(0.3, 0.3, 100.0, 1.0e18, {0.0, 0.0}); // electron density in #/m^3
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
	delete paravis;
	delete paramesh;

	// this is a testing section
		// solve the poisson equation using a finite difference method
	//Simulation mysim = Simulation(FINITE_DIFFERENCE, paramesh);
	
	/*
	// loop over the internal nodes adding the laplace operator to the matrix
	double ** mat = new double *[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x(); i++) mat[i] = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=0; i<paramesh->reg_num_nodes_x(); i++){ // columns
		for (auto j=0; j<paramesh->reg_num_nodes_y(); j++) mat[i][j] = 0;
	}

	unsigned int cind, lind, rind, uind, dind;
	for (auto i=1; i<paramesh->reg_num_nodes_x()-1; i++){ // columns
		for (auto j=1; j<paramesh->reg_num_nodes_y()-1; j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			lind = paramesh->reg_inds_to_glob_ind(i-1,j);
			rind = paramesh->reg_inds_to_glob_ind(i+1,j);
			uind = paramesh->reg_inds_to_glob_ind(i,j+1);
			dind = paramesh->reg_inds_to_glob_ind(i,j-1);
			mat[cind][cind] = -4;
			mat[cind][lind] = 1;
			mat[cind][rind] = 1;
			mat[cind][uind] = 1;
			mat[cind][dind] = 1;
		}
	}

	// write in the boundary conditions
	// zeros for now...just testing the code

	// construct the right side using the density field
	double q_electron = -1.6e-19, eps0 = 8.854e-12;
	const double * reldata = &paramesh->data("e_density");
	double * rhs = new double[paramesh->reg_num_nodes_y()*paramesh->reg_num_nodes_x()];
	for (auto i=1; i<paramesh->reg_num_nodes_x()-1; i++){ // columns
		for (auto j=1; j<paramesh->reg_num_nodes_y()-1; j++){ // rows
			cind = paramesh->reg_inds_to_glob_ind(i,j);
			rhs[cind] = -reldata[cind]*q_electron/eps0;
		}
	}
	*/

	// now solve the system (PETSc)

	// visualize the solution by putting it in the mesh




	/* this is what I would do if I had all the classes ready
	// create an equation
	// del^2(potential) = -e*(e_density)/epsilon
	//Equation poiss = Equation(EQUATION_POISSON);

	// set boundary conditions


	// create a simulation

	// run the simulation

	// visualize the results
	*/

	return 0;
}