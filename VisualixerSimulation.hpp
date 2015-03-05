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

	// inspectors


	// mutators
	void bind_simulation(const SimulationData & simdat);
	void set_frequency_Hz(unsigned int freq); // set the # of cycles per second
	void set_colorby_field(std::string fieldname);
	void set_alpha_field(std::string fieldname);


	// action
	void run();

protected:

private:

	void increment_time_step();

	void onPrepareData();
	void onRender();
	bool MainLoop();
	void onExit();

	GLuint * line_elements;
	GLuint lebo;

	const SimulationData * _simdata;

	std::string _colorby_field, _alpha_field;
	unsigned int _freq_Hz, _cur_time_step; // frequency at which the simulation plays
	unsigned int num_line_elements, num_per_line_element, line_element_offset;
};

#endif
