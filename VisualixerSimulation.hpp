#ifndef _VISUALIXERSIMULATION_H
#define _VISUALIXERSIMULATION_H

#include "SimulationData.hpp"
#include "Visualixer.hpp"


#include <iostream>
#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>

// this needs to be able to view static images or "movies" from a simulation
// One of the options for a static image in 3D could be slices in x, y, z
// ... this could also be a dynamic image option
//*********** Simulation Data Visualixer ****************
class simulation_visualixer : public visualixer{
public:
	simulation_visualixer();
	~simulation_visualixer();

	void add_simulation(const SimulationData & simdat);

protected:
	void onRender();

	bool MainLoop();
	void onExit();

private:
	const SimulationData * _simdata;

};

#endif
