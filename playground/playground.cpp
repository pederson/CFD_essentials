#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"
#include "../FDTDSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double dx = 5.0e-4, srcfreq = 15.0e+13;
	double srclocx = 1.5e-6, srclocy = 5.0e-6, srclocz = 1.5e-6;
	unsigned int numiters = 300;
	//double dt = 0.5*dx/c0;

	// ring model
	/*
	srcfreq = 1.0e+9;
	srclocx = 0.1; 
	srclocy = 0.5;
	ParametricModel2D paramodel;
	paramodel.set_model_name("Ring");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	Circle c1 = Circle(0.1, {0.5, 0.5}, paramodel.material("Dielectric"));
	Circle c2 = Circle(0.06, {0.5, 0.5}, paramodel.material("Vacuum"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.005;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/
	

	// holey waveguide model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("Holey Waveguide");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	rectangle r1 = rectangle(9.0e-6, 2.0e-6, {4.5e-6, 3.0e-6}, paramodel.material("Dielectric"));
	circle c1 = circle(0.5e-6, {1.5e-6, 3.0e-6}, paramodel.material("Vacuum"));
	circle c2 = circle(0.5e-6, {3.5e-6, 3.0e-6}, paramodel.material("Vacuum"));
	circle c3 = circle(0.5e-6, {5.5e-6, 3.0e-6}, paramodel.material("Vacuum"));
	circle c4 = circle(0.5e-6, {7.5e-6, 3.0e-6}, paramodel.material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	//*/


	// holey waveguide defect model
	/*
	ParametricModel2D paramodel;
	paramodel.set_model_name("Holey Defect");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	Rectangle r1 = Rectangle(1.0, 0.2, {0.5, 0.5}, paramodel.material("Dielectric"));
	Circle c1 = Circle(0.05, {0.11, 0.5}, paramodel.material("Vacuum"));
	Circle c2 = Circle(0.05, {0.31, 0.5}, paramodel.material("Vacuum"));
	Circle c3 = Circle(0.05, {0.59, 0.5}, paramodel.material("Vacuum"));
	Circle c4 = Circle(0.05, {0.79, 0.5}, paramodel.material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	paramodel.add_object(&c3);
	paramodel.add_object(&c4);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.005;
	srcfreq = 1.0e+9;
	srclocx = 0.03;
	srclocy = 0.5;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.2, 0.8, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// rod array model
	/*
	ParametricModel2D paramodel;
	paramodel.set_model_name("Rod Array");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {6.0});
	Circle c1 = Circle(0.5e-6, {9.5e-6, 1.0e-6}, paramodel.material("Dielectric"));
	paramodel.create_lattice(&c1, {2.0e-6, 0}, {0, 2.0e-6}, 5, 5);
	paramodel.print_summary();
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.5e-7;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 3.0e-5, 0.0, 1.0e-5, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// current rod model
	/*
	parametric_model_2d paramodel;
	paramodel.set_model_name("Rod With Current");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_physical_property("current_density");
	paramodel.add_material("Vacuum", {1.0, 0.0});
	paramodel.add_material("Dielectric", {6.0, 1.0e+8});
	paramodel.add_material("DielectricInverse", {6.0, -1.0e+8});
	circle c1 = circle(0.005, {0.03, 0.05}, paramodel.material("Dielectric"));
	circle c2 = circle(0.005, {0.07, 0.05}, paramodel.material("DielectricInverse"));
	paramodel.add_object(&c1);
	paramodel.add_object(&c2);
	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 0.1, 0.0, 0.1, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/


	// metal aperture model
	/*
	srcfreq = 1.0e+10;
	srclocx = 0.01; 
	srclocy = 0.05;
	ParametricModel2D paramodel;
	paramodel.set_model_name("Aperture");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_physical_property("current_density");
	paramodel.add_material("Vacuum", {1.0, 0.0});
	paramodel.add_material("Dielectric", {1000.0, 0.0});
	Rectangle r1 = Rectangle(0.3/40, 0.2, {0.05, 0.05}, paramodel.material("Dielectric"));
	Rectangle r2 = Rectangle(0.3/40, 0.03, {0.05, 0.05}, paramodel.material("Vacuum"));
	paramodel.add_object(&r1);
	paramodel.add_object(&r2);
	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 0.1, 0.0, 0.1, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// conductive slab model
	/*
	srclocx = 1.5e-6, srclocy = 5.0e-6;
	ParametricModel2D paramodel;
	paramodel.set_model_name("ConductiveSlab");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_physical_property("conductivity");
	paramodel.add_material("Vacuum", {1.0, 0.0});
	paramodel.add_material("Dielectric", {1.5, 10000.5});
	Rectangle r1 = Rectangle(1.5e-5, 1.0e-5, {1.0e-5, 5.0e-6}, paramodel.material("Dielectric"));
	paramodel.add_object(&r1);
	paramodel.print_summary();
	// convert the model into a mesh
	dx = 0.5e-7;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 2.0e-5, 0.0, 1.0e-5, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// frequency dependent slab model
	/*
	dx = 1.0e-6;
	srcfreq = 3.0e+12;	// should be reflected
	//srcfreq = 3.0e+13;	// should pass right through
	srclocx = 10*dx;
	srclocy = 50*dx;
	ParametricModel2D paramodel;
	paramodel.set_model_name("FrequencyDepenedentSlab");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_physical_property("conductivity");
	paramodel.add_physical_property("polenumerator");
	paramodel.add_physical_property("polefreq");
	paramodel.add_material("Vacuum", {1.0, 0.0, 0.0, 0.0});
	paramodel.add_material("Dielectric", {1.0, 349.5e+3, -39.47e+15, 1.0e+11});
	Rectangle r1 = Rectangle(200*dx, 100*dx, {200*dx, 50*dx}, paramodel.material("Dielectric"));
	paramodel.add_object(&r1);
	paramodel.print_summary();
	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 500*dx, 0.0, 100*dx, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// 1D free space
	/*
	dx = 0.005;
	srcfreq = 1.0e+9;
	srclocx = 0.3;
	RegularMesh paramesh;
	paramesh = RegularMesh::create_regular_grid_b(dx, 0.0, 1.0);
	//*/

	// 2D free space model
	/*
	srcfreq = 1.0e+9;
	srclocx = 0.5; 
	srclocy = 0.5;
	ParametricModel2D paramodel;
	paramodel.set_model_name("EmptySpace");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {1.0});
	Circle c1 = Circle(0.1, {0.5, 0.5}, paramodel.material("Vacuum"));
	paramodel.add_object(&c1);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.005;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// 2D lens
	
	srcfreq = 6.0e+9;
	srclocx = 0.1; 
	srclocy = 0.5;
	ParametricModel2D paramodel;
	paramodel.set_model_name("Lens2D");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {3.0});
	Parabola c1 = Parabola({0.3, 0.5}, {0.5, 0.5}, 0.2, paramodel.material("Dielectric"));
	paramodel.add_object(&c1);
	Parabola c2 = Parabola({0.7, 0.5}, {0.5, 0.5}, 0.2, paramodel.material("Dielectric"));
	paramodel.add_object(&c2);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.005;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	numiters = 800;
	//*/

	// 3D free space
	/*
	dx = 0.02;
	srcfreq = 1.0e+9;
	srclocx = 0.4;
	srclocy = 0.4;
	srclocz = 0.3;
	RegularMesh paramesh;
	paramesh = RegularMesh::create_regular_grid_b(dx, 0.0, 1.0, 0.0, 0.8, 0.0, 0.6);
	//*/

	// 3D dielectric sphere
	/*
	srcfreq = 1.0e+9;
	srclocx = 0.2;
	srclocy = 0.2;
	srclocz = 0.2;
	ParametricModel3D paramodel;
	paramodel.set_model_name("DielectricSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {3.0});
	Sphere s1 = Sphere(0.2, {0.5, 0.5, 0.5}, paramodel.material("Dielectric"));
	paramodel.add_object(&s1);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.015;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_3d(paramodel, dx, 0.0, 1.0, 0.0, 0.8, 0.0, 0.6, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// 3D dielectric cylinder
	/*
	srcfreq = 1.0e+9;
	srclocx = 0.2;
	srclocy = 0.2;
	srclocz = 0.2;
	ParametricModel3D paramodel;
	paramodel.set_model_name("DielectricSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {3.0});
	Cylinder s1 = Cylinder(0.2, 0.5, {1.0, 1.0, 0.0}, {0.5, 0.5, 0.5}, paramodel.material("Dielectric"));
	paramodel.add_object(&s1);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.015;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_3d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// 3D dielectric box
	/*
	srcfreq = 1.0e+9;
	srclocx = 0.2;
	srclocy = 0.2;
	srclocz = 0.2;
	ParametricModel3D paramodel;
	paramodel.set_model_name("DielectricSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {3.0});
	Box s1 = Box(0.5, 0.2, 0.3, {0.0, 0.0, 1.0}, {0.5, 0.5, 0.5}, paramodel.material("Dielectric"));
	paramodel.add_object(&s1);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.015;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_3d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/

	// 3D dielectric lens
	/*
	srcfreq = 5.0e+9;
	srclocx = 0.1;
	srclocy = 0.5;
	srclocz = 0.5;
	ParametricModel3D paramodel;
	paramodel.set_model_name("DielectricSphere");
	paramodel.add_physical_property("eps_rel");
	paramodel.add_material("Vacuum", {1.0});
	paramodel.add_material("Dielectric", {3.0});
	ParabolicDish s1 = ParabolicDish({0.35, 0.5, 0.5}, {0.55, 0.5, 0.5}, 0.15, paramodel.material("Dielectric"));
	paramodel.add_object(&s1);
	ParabolicDish s2 = ParabolicDish({0.65, 0.5, 0.5}, {0.45, 0.5, 0.5}, 0.15, paramodel.material("Dielectric"));
	paramodel.add_object(&s2);
	// convert the model into a mesh
	cout << "about to make mesh" << endl;
	dx = 0.015;
	RegularMesh paramesh;
	paramesh = build_simple_mesh_3d(paramodel, dx, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, paramodel.material("Vacuum"));
	paramesh.print_summary();
	//*/


	
	
	// view the mesh
	mesh_visualixer paravis;
	paravis.bind_mesh(paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	paravis.set_colorby(&paramesh.data("eps_rel"));
	paravis.run();
	//*/
	


	// initialize the simulation
	FDTDSimulation fsim;
	fsim.bind_mesh(paramesh);
	fsim.bind_rel_permittivity(&paramesh.data("eps_rel"));
	//fsim.bind_conductivity(&paramesh.data("conductivity"));
	//fsim.bind_single_pole(&paramesh.data("polenumerator"), &paramesh.data("polefreq"));
	//fsim.bind_current_density_z(&paramesh.data("current_density"));
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy + 0.1, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy - 0.1, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy + 0.2, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy - 0.2, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy + 0.3, srclocz);
	fsim.add_sinusoidal_source(srcfreq, 0.0, srclocx, srclocy - 0.3, srclocz);

	//fsim.add_sinusoidal_source(15.0e+13, 0.0, 1.5e-6, 2.5e-6);
	//fsim.add_gaussian_source(10.0, 10.0, 4.5e-6, 4.0e-6);
	fsim.set_num_iters(numiters);
	fsim.run();
	fsim.view_results();
	fsim.output_HDF5("plasmaslab.h5");

	/*
	SimulationData fsimdat = SimulationData::read_HDF5("../../testfiles/plasmaslab.h5");
	simulation_visualixer vsim;
	vsim.bind_simulation(fsimdat);
	vsim.set_colorby_field("E_z");
	vsim.run();
	*/

	return 0;
}
