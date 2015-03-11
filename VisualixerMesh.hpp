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

	void bind_mesh(const Mesh & mesh);
	void set_test_case();

protected:

	void onPrepareData();

	const Mesh * _mesh;
};

#endif