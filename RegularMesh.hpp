#ifndef _REGULARMESH_H
#define _REGULARMESH_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <fstream>

#include "Mesh.hpp"


// regular mesh
class RegularMesh : public Mesh {
public:

  //RegularMesh();
  //RegularMesh(const RegularMesh & rmesh);
  //~RegularMesh();

  // specialty stuff for regular grids
  MeshNode & regular_node(unsigned int i, unsigned int j=0, unsigned int k=0);
  unsigned int reg_num_nodes_x() const {return _num_nodes_x;};
  unsigned int reg_num_nodes_y() const {return _num_nodes_y;};
  unsigned int reg_num_nodes_z() const {return _num_nodes_z;};
  //unsigned int & reg_nodes_boundary(BoundaryLocation loc);
  std::vector<unsigned int> reg_left_inds() const;
  std::vector<unsigned int> reg_right_inds() const;
  std::vector<unsigned int> reg_top_inds() const;
  std::vector<unsigned int> reg_bottom_inds() const;
  std::vector<unsigned int> reg_front_inds() const;
  std::vector<unsigned int> reg_back_inds() const;
  unsigned int reg_inds_to_glob_ind(unsigned int i, unsigned int j=0, unsigned int k=0) const;
  unsigned int nearest_node(double x_loc, double y_loc=0.0, double z_loc=0.0) const;

  // grid generation and refinement
  static RegularMesh create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static RegularMesh create_regular_grid_b(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);

private:
  unsigned int _num_nodes_x, _num_nodes_y, _num_nodes_z;
  double _res;

  void create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen=0.0, double ycen=0.0, double zcen=0.0);


};



#endif
