#ifndef _VISUALIXERMESH_H
#define _VISUALIXERMESH_H

#include "Mesh.hpp"
#include "Visualixer.hpp"


#include <iostream>
#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>

class mesh_visualixer : public visualixer{
public:
	mesh_visualixer();
	~mesh_visualixer();

	void add_mesh(Mesh * mesh);
	void set_test_case();

protected:
	void onRender();

	bool MainLoop();
	void onExit();

	GLuint * line_elements;
	GLuint lebo;
	unsigned int num_line_elements, num_per_line_element, line_element_offset;
};

#endif