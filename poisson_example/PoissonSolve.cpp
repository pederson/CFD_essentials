#include "../Visualixer.hpp" 
#include "../GeometricObject.hpp"
#include "../Mesh.hpp" 
#include "../Converter.hpp" 
//#include "OperatorTerm.hpp"
//#include "Equation.hpp"
//#include "Simulation.hpp"

int main(int argc, char * argv[]){

	// define the geometry model that will be used
	parametric_model_2d my_param2;
	Static_Mesh * paramesh;
	my_param2.set_model_name("GaussianDistElectrons");
	my_param2.add_physical_property("e_density");
	my_param2.add_material("Air", {1.0});
	my_param2.add_material("Dielectric", {5.0});
	my_param2.add_material("Dielectric2", {9.0});
	gaussian_2d ga1 = gaussian_2d(0.3, 0.3, 10.0, 1.0, {0.0, 0.0});
	my_param2.add_object(&ga1);

	// convert the model into a mesh
	paramesh = build_simple_mesh_2d(&my_param2, 0.01, -1.0, 1.0, -1.0, 1.0, my_param2.get_material("Air"));
	

	// view the mesh
	mesh_visualixer * paravis = new mesh_visualixer();
	paravis->add_mesh(paramesh);
	paravis->set_color_ramp(CRamp::DIVERGENT_9);
	paravis->set_colorby(&paramesh->data("e_density"));
	paravis->run();
	delete paravis;
	delete paramesh;

	// create an equation
	// -del^2(potential) = e*(e_density)/epsilon

	// set boundary conditions

	// start a simulation

	// run the simulation

	// visualize the results

	return 0;
}