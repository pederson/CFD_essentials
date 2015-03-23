#include "../VisualixerMesh.hpp"
#include "../RegularMesh.hpp" 
#include "../VisualixerSimulation.hpp"
#include "../Converter.hpp"
#include "../PlasmaDriftDiffusionSimulation.hpp"
#include "../PlasmaMomentumSimulation.hpp"
#include "../FDTDSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){
	// constants
	double dx = 0.5*1.36e-5;
	//double dt = 1.0e16*9.0e-12;
	double plasma_dens = 1.0e+22;		// in #/m^3
	plasma_dens = 2;

	// dot of plasma model
	ParametricModel2D paramodel;
	paramodel.set_model_name("PlasmaDot");
	paramodel.add_physical_property("density");
	paramodel.add_physical_property("init_vel");
	paramodel.add_material("Vacuum", {0.0, 0.0});
	paramodel.add_material("Plasma", {plasma_dens, 0.0});
	Circle c1 = Circle(6.0e-5, {2.7e-3/2, 2.7e-3/2}, paramodel.get_material("Plasma"));
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


	// initialize the simulations
	PlasmaMomentumSimulation msim;
	PlasmaDriftDiffusionSimulation dsim;
	FDTDSimulation fdtdsim;
	unsigned int numiters = 300;

	// Diffusion simulation
	dsim.bind_mesh(paramesh);
	dsim.set_initial_density(&paramesh.data("density"));
	//dsim.set_time_step(dt);
	dsim.set_num_iters(numiters);

	// Momentum Simulation
	msim.bind_mesh(paramesh);
	msim.set_num_iters(numiters);
	msim.set_initial_velocity_x(&paramesh.data("init_vel"));
	msim.set_initial_velocity_y(&paramesh.data("init_vel"));
	msim.bind_E_z(fdtdsim.E_z_ptr());
	msim.bind_collision_rate(dsim.density_ptr());

	// FDTD Simulation
	cvector eps_rel;
	eps_rel.multiply(dsim.density_ptr());
	eps_rel.set_additive_constant(1.0);
	if (eps_rel.isempty()) cout << "eps rel is empty" << endl;
	fdtdsim.bind_mesh(paramesh);
	fdtdsim.bind_rel_permittivity(eps_rel);
	//fdtdsim.bind_current_density_z(msim.velocity_z_ptr());
	fdtdsim.add_sinusoidal_source(3.0e+12/6.0, 0.0, 1.0e-3, 1.0e-3);
	fdtdsim.set_num_iters(numiters);
	fdtdsim.run();

	for (auto i=0; i<numiters; i++){
		//dsim.run(1);
		//msim.run(1);
		//fdtdsim.run(1);
	}

	fdtdsim.view_results();




	return 0;
}
