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
  unsigned int core_group;

  std::vector<unsigned int> neighbor_keys;		// keys of neighbor points

  // container for physical properties
  std::vector<double> phys_properties;

  void print_summary();

protected:

private:

};

// cell class for finite volume methods
class Mesh_Cell{
public:
  unsigned int num_vertices;
  std::vector<unsigned int> vertices;

protected:

private:

};

class Mesh_Element{

};

class Static_Mesh{
public:
  Static_Mesh();   // constructor
  ~Static_Mesh();    // destructor

  void print_summary() const;
  void print_detailed() const;
  

  // metadata access
  MeshType get_mesh_type() const {return _mesh_type;};
  unsigned int get_num_dims() const {return _num_dims;};
  unsigned int nodecount() const {return 0;};
  unsigned int elementcount() const {return 0;};

  double x_min() const {return _xmin;};
  double y_min() const {return _ymin;};
  double z_min() const {return _zmin;};
  double x_max() const {return _xmax;};
  double y_max() const {return _ymax;};
  double z_max() const {return _zmax;};


  // node access
  const Node * node_ptr(unsigned int i) const;

  // property interaction
  void add_phys_property(std::string property_name, double * property_vals);
  
  void set_background_properties(std::vector<double> properties);
  float * get_phys_property_ptr(std::string property_name);


  // grid generation and refinement
  static Static_Mesh * create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static Static_Mesh * create_regular_grid(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);
  //static Mesh * create_unstructured_tri_simple();

protected:

private:
  // metadata
  MeshType _mesh_type;
  unsigned int _num_dims;
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  // nodes and keys
  Node * _nodes; // contains keys to the nodes

  // user-defined propertie for the mesh
  std::vector<std::string> _phys_property_names; // the name position in this vector corresponds with the position of the property value in the node
  std::map<std::string, double *> _properties;

  void calc_extents();

};

class Mutable_Mesh{
public:

  Mutable_Mesh();		// constructor
  ~Mutable_Mesh();		// destructor

  void print_summary();
  

  // member data accessors
  MeshType get_mesh_type() const {return mesh_type;};
  unsigned int get_num_dims() const {return num_dims;};
  unsigned int get_num_nodes() const {return mesh_nodes.size();};
  void set_mesh_type(MeshType type);
  void set_num_dims(unsigned int ndims);
  void set_num_nodes(unsigned int number_of_nodes);

  double get_xmin(){return xmin;};
  double get_ymin(){return ymin;};
  double get_zmin(){return zmin;};
  double get_xmax(){return xmax;};
  double get_ymax(){return ymax;};
  double get_zmax(){return zmax;};
  void set_xmin(double x_min);
  void set_ymin(double y_min);
  void set_zmin(double z_min);
  void set_xmax(double x_max);
  void set_ymax(double y_max);
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
  void set_background_properties(std::vector<double> properties);
  float * get_phys_property_ptr(std::string property_name);


  // grid generation and refinement
  static Mutable_Mesh * create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static Mutable_Mesh * create_regular_grid(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
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

  void calc_extents();
};




// unstructured cell transformation here


#endif
