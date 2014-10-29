#ifndef _MESH_CLASS_H
#define _MESH_CLASS_H

#include <stdlib>
#include <vector>
#include <list>

#include "node_class.hpp"

enum MeshType{REGULAR=0, UNSTRUCTURED_TRI=1, UNSTRUCTURED_QUAD=2};

struct Node{
  unsigned int index; 		// index of point within x,y,z data
  unsigned int num_neighbors;	// number of neighbors
  std::vector<unsigned int> neighbor_index;		// index of neighbor point within x,y,z data

  // physical properties
  double epsilon; // relative dielectric constant
  double mu; // relative permissivity constant
};

// cell class for finite volume methods
class Cell{
public:
  unsigned int num_vertices;
  std::vector<unsigned int> vertices;

protected:

private:

};

class Mesh{
public:

  Mesh();		// constructor
  ~Mesh();		// destructor

  create_regular_grid(unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z); // create a regular grid of points and store it in the mesh


protected:

private:
  MeshType mesh_type;
  unsigned int num_nodes; // number of nodes in the mesh
  double *x, *y, *z; // the coordinates of each node (z optional)
  bool *boundary; // flag for boundary points

  unsigned int *core_group; // which core does this node belong to (for parallel processing)

  std::list<Node> mesh_nodes; // contains Node structures
};


// unstructured cell transformation here


#endif
