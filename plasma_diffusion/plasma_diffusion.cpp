#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../SimulationData.hpp"
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"
#include "../PlasmaDriftDiffusionSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double dx = 0.5*1.36e-5;
	//double dt = 1.0e16*9.0e-12;
	double plasma_dens = 1.0e+22;		// in #/m^3

	// dot of plasma model
	
	parametric_model_2d paramodel;
	paramodel.set_model_name("PlasmaDot");
	paramodel.add_physical_property("density");
	paramodel.add_material("Vacuum", {0.0});
	paramodel.add_material("Plasma", {plasma_dens});
	circle c1 = circle(6.0e-5, {2.7e-3/2, 2.7e-3/2}, paramodel.get_material("Plasma"));
	paramodel.add_object(&c1);
	//*/

	// convert the model into a mesh
	RegularMesh paramesh;
	paramesh = build_simple_mesh_2d(paramodel, dx, 0.0, 2.7e-3, 0.0, 2.7e-3, paramodel.get_material("Vacuum"));
	paramesh.print_summary();

	// view the mesh
	mesh_visualixer paravis;
	paravis.bind_mesh(paramesh);
	paravis.set_color_ramp(CRamp::DIVERGENT_9);
	paravis.set_colorby(&paramesh.data("density"));
	paravis.run();


	// initialize the simulation
	PlasmaDriftDiffusionSimulation psim;
	psim.bind_mesh(paramesh);
	psim.set_initial_density(&paramesh.data("density"));
	//psim.set_time_step(dt);
	psim.set_num_iters(800);
	psim.run();
	psim.view_results();
	//psim.output_HDF5("testout.h5");


	return 0;
}