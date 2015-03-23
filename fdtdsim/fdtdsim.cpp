#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"
#include "../FDTDSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double dx = 0.03e-6;
	//double dt = 0.5*dx/c0;

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
	paramodel.set_model_name("Ring");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	circle c1 = circle(0.1, {0.5, 0.5}, paramodel.get_material("Dielectric"));
	circle c2 = circle(0.06, {0.5, 0.5}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	//*/
	

	// holey waveguide model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("Holey Waveguide");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	rectangle r1 = rectangle(9.0e-6, 2.0e-6, {4.5e-6, 3.0e-6}, paramodel.get_material("Dielectric"));
	circle c1 = circle(0.5e-6, {1.5e-6, 3.0e-6}, paramodel.get_material("Vacuum"));
	circle c2 = circle(0.5e-6, {3.5e-6, 3.0e-6}, paramodel.get_material("Vacuum"));
	circle c3 = circle(0.5e-6, {5.5e-6, 3.0e-6}, paramodel.get_material("Vacuum"));
	circle c4 = circle(0.5e-6, {7.5e-6, 3.0e-6}, paramodel.get_material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	//*/

	// holey waveguide defect model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("Holey Defect");
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

	// rod array model

	ParametricModel2D paramodel;
	paramodel.set_model_name("Rod Array");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_physical_property("current_density_z");
	paramodel.add_material("Vacuum", {1.0, 0.0});
	paramodel.add_material("Dielectric", {6.0, 0.0});
	paramodel.add_material("CurrentWire", {10.0, 1.0e+7});
	Circle c1 = Circle(0.1e-6, {1.5e-6, 3.0e-6}, paramodel.get_material("CurrentWire"));
	Circle c2 = Circle(0.5e-6, {3.5e-6, 4.0e-6}, paramodel.get_material("Dielectric"));
	Circle c3 = Circle(0.5e-6, {5.5e-6, 4.0e-6}, paramodel.get_material("Dielectric"));
	Circle c4 = Circle(0.5e-6, {7.5e-6, 4.0e-6}, paramodel.get_material("Dielectric"));
	//circle c5 = circle(0.5e-6, {1.5e-6, 2.0e-6}, paramodel.get_material("Dielectric"));
	Circle c6 = Circle(0.5e-6, {3.5e-6, 2.0e-6}, paramodel.get_material("Dielectric"));
	Circle c7 = Circle(0.5e-6, {5.5e-6, 2.0e-6}, paramodel.get_material("Dielectric"));
	Circle c8 = Circle(0.5e-6, {7.5e-6, 2.0e-6}, paramodel.get_material("Dielectric"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	//paramodel.add_object(&c5);
	paramodel.add_object(&c6);
	paramodel.add_object(&c7);
	paramodel.add_object(&c8);
	//*/
	

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 9.0e-6, 0.0, 6.0e-6, paramodel.get_material("Vacuum"));
	paramesh.print_summary();

	// view the mesh
	mesh_visualixer paravis;
	paravis.bind_mesh(paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	paravis.set_colorby(&paramesh.data("eps_rel"));
	paravis.run();


	// initialize the simulation
	FDTDSimulation fsim;
	fsim.bind_current_density_z(&paramesh.data("current_density_z"));
	fsim.bind_mesh(paramesh);
	fsim.bind_rel_permittivity(&paramesh.data("eps_rel"));
	fsim.add_sinusoidal_source(3.0e+14/6.0, 0.0, 0.5e-6, 3.0e-6);
	//fsim.add_gaussian_source(10.0, 10.0, 4.5e-6, 4.0e-6);
	fsim.set_num_iters(800);
	fsim.run();
	fsim.view_results();
	//fsim.output_HDF5("testout.h5");


	return 0;
}
