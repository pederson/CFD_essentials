#ifndef _MESH_CLASS_H
#define _MESH_CLASS_H

#include <stdlib.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>

enum MeshType{REGULAR=0, UNSTRUCTURED_TRI=1, UNSTRUCTURED_QUAD=2};

class Node{
public:

  Node();
  ~Node();

  unsigned int key; // key associated with this node

  double x, y, z; 
  bool boundary;

  std::vector<unsigned int> neighbor_keys;		// keys of neighbor points

  unsigned int core_group;

  // container for physical properties
  std::vector<double> phys_properties;

  void print_summary();

protected:

private:

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

  void print_summary();
  void calc_extents();

  // member data accessors
  MeshType get_mesh_type(){return mesh_type;};
  void set_mesh_type(MeshType type);
  unsigned int get_num_dims(){return num_dims;};
  void set_num_dims(unsigned int ndims);
  unsigned int get_num_nodes(){return mesh_nodes.size();};
  void set_num_nodes(unsigned int number_of_nodes);

  double get_xmin(){return xmin;};
  void set_xmin(double x_min);
  double get_ymin(){return ymin;};
  void set_ymin(double y_min);
  double get_zmin(){return zmin;};
  void set_zmin(double z_min);
  double get_xmax(){return xmax;};
  void set_xmax(double x_max);
  double get_ymax(){return ymax;};
  void set_ymax(double y_max);
  double get_zmax(){return zmax;};
  void set_zmax(double z_max);

  // node access and manipulation
  Node * get_node_ptr(unsigned int key);
  unsigned int get_node_key(unsigned int i);
  //unsigned int get_node_index(unsigned int key);
  void add_node(Node * new_node); // add node and add neighbor connections
  void remove_node(unsigned int key); // remove node and delete neighbor connections

  // property interaction
  void add_phys_property(std::string property_name);
  //void set_phys_property(std::string property_name, double property_value);
  unsigned int get_phys_property_position(std::string property_name);


  // grid generation and refinement
  static Mesh * create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static Mesh * create_regular_grid(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);
  //static Mesh * create_unstructured_tri_simple();

protected:

private:
  // auxillary information
  MeshType mesh_type;
  unsigned int num_dims;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  // nodes and keys
  std::vector<unsigned int> node_keys; // contains keys to the nodes
  std::map<unsigned int, Node *> mesh_nodes; // contains pointers to Node structures as a list (so that the mesh is refinable)

  // physical properties on the mesh
  std::vector<std::string> phys_property_names; // the name position in this vector corresponds with the position of the property value in the node
};


// unstructured cell transformation here


#endif
