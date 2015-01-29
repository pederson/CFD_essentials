#include "../VisualixerMesh.hpp" 
#include "../GeometricObject.hpp"
#include "../Mesh.hpp" 
#include "../Converter.hpp" 

int main(int argc, char * argv[]){
	// test a mesh viewer made from a parametric model
	mesh_visualixer * paravis = new mesh_visualixer();
	parametric_model_2d my_param2;
	Static_Mesh * paramesh;
	my_param2.set_model_name("Display Test Circle");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});
	my_param2.add_material("Dielectric2", {9.0, 2.0});
	circle circout = circle(0.4, {0.0, 0.0}, my_param2.get_material("Dielectric"));
	circle circin = circle(0.3, {0.0, 0.0}, my_param2.get_material("Air"));
	rectangle gap = rectangle(0.05, 0.1, {0.0, 0.35}, my_param2.get_material("Air"));
	my_param2.add_object(&circout);
	my_param2.add_object(&circin);
	my_param2.add_object(&gap);
	paramesh = build_simple_mesh_2d(&my_param2, 0.02, -1.0, 1.0, -1.0, 1.0, my_param2.get_material("Air"));
	paravis->add_mesh(paramesh);
	paravis->set_color_ramp(CRamp::DIVERGENT_9);
	paravis->set_colorby(&paramesh->data("Epsilon_rel"));
	//paravis->set_colorby(&paramesh->num_connections());
	paravis->run();
	delete paravis;
	delete paramesh;

	return 0;
}