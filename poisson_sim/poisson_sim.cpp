//#include "../Visualixer.hpp" 
#include "../VisualixerMesh.hpp"
#include "../GeometricObject.hpp"
#include "../RegularMesh.hpp" 
#include "../Converter.hpp" 
#include "../PoissonSimulation.hpp"


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
	double dx = 0.002;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(my_param2, dx, -0.1, 0.1, -0.15, 0.15, my_param2.get_material("Air"));
	
	// view the mesh
	mesh_visualixer paravis;
	paravis.bind_mesh(paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	paravis.set_colorby(&paramesh.data("e_density"));
	paravis.run();
	//delete paramesh; // WHY DOES THIS COMPILE IF THIS IS UNCOMMENTED???


	// construct the right side using the density field
	double q_electron = -1.6e-19, eps0 = 8.854e-12, m_electron = 9.11e-31;

	cvector rightside;
	rightside.multiply(&paramesh.data("e_density"));
	rightside.multiply(q_electron);
	rightside.divide(eps0*dx*dx);

	PoissonSimulation poissim;
	poissim.bind_mesh(paramesh);
	poissim.bind_rhs(rightside);
	poissim.run();
	poissim.view_results();

	
	return 0;
}