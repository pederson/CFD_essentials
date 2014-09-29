#ifndef _MESH_CLASS_H
#define _MESH_CLASS_H

#include <stdlib>
#include "node_class.hpp"

// MESH class definitions go here 
class Mesh{
	public:

		unsigned int meshtype;
		unsigned int num_nodes;
		std::list(Node) nodes;
		std::list(std::vector<int>) neighbors;
}


#endif