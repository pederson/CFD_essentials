#include "../../RegularMesh.hpp" 
#include "../../SimulationData.hpp"
#include "../../VisualixerSimulation.hpp"

using namespace std;

int main(int argc, char * argv[]){

	
	SimulationData fsimdat = SimulationData::read_HDF5(argv[1]);
	vector<string> fieldsnames = fsimdat.fieldnames();

	simulation_visualixer vsim;
	vsim.bind_simulation(fsimdat);
	vsim.set_colorby_field(fieldsnames.at(0));
	vsim.set_color_ramp(CRamp::MATLAB_PARULA);
	vsim.set_snapshot_increment(3);
	vsim.run();
	

	return 0;
}
