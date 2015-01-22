#ifndef _MESH_H
#define _MESH_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <fstream>


enum MeshType{REGULAR=0, UNSTRUCTURED_TRI=1, UNSTRUCTURED_QUAD=2};
enum ElementType{EMPTY, POINT, LINE, TRIANGLE, QUADRANGLE, TETRAHEDRON, UNKNOWN};

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


class Mesh_Node{
public:
  Mesh_Node(double x, double y, double z=0.0, bool boundary=false, unsigned int num_connections=0, unsigned int core_group=false);
  ~Mesh_Node();

  // utils
  void print_summary() const;
  void print_detailed() const;
  static double dist_sq(Mesh_Node * node1, Mesh_Node * node2);
  static double area_tri(Mesh_Node * node1, Mesh_Node * node2, Mesh_Node * node3);

  // member data access
  double x() const {return _x;};
  double y() const {return _y;};
  double z() const {return _z;};
  unsigned int num_connections() const {return _num_connections;};
  unsigned int core_group() const {return _core_group;};
  bool boundary() const {return _boundary;};

  // mutators
  void add_connection() {_num_connections++;};
  void remove_connection() {_num_connections--;};
  void set_x(double xval) {_x = xval;};
  void set_y(double yval) {_y = yval;};
  void set_z(double zval) {_z = zval;};
  void set_boundary(bool boundary) {_boundary = boundary;};
  void set_core_group(unsigned int core_group) {_core_group = core_group;};
  void set_num_connections(unsigned int num_connections) {_num_connections = num_connections;};

private:
  double _x, _y, _z; 
  bool _boundary;
  unsigned int _core_group, _num_connections;

};


class Mesh_Element{
public:
  Mesh_Element(std::vector<unsigned int> vertex_inds);
  ~Mesh_Element();

  // utils
  void print_summary() const;
  void print_detailed() const;
  double area() const;
  double perimeter() const;

  // member data access
  unsigned int num_vertices() const {return _vertex_inds.size();};
  std::vector<unsigned int> vertex_inds() const {return _vertex_inds;};

  // mutators
  void remove_vertex(unsigned int vert_ind);
  void add_vertex(unsigned int vert_ind, int position=-1);

private:
  std::vector<unsigned int> _vertex_inds;
  ElementType _element_type; // 

  void recalc_type();

};


class Static_Mesh{
public:
  Static_Mesh();   // constructor
  ~Static_Mesh();    // destructor

  void print_summary() const;
  void print_detailed() const;
  

  // metadata access
  MeshType mesh_type() const {return _mesh_type;};
  unsigned int num_dims() const {return _num_dims;};
  unsigned int nodecount() const {return _nodecount;};
  unsigned int elementcount() const {return _elementcount;};

  double xmin() const {return _xmin;};
  double ymin() const {return _ymin;};
  double zmin() const {return _zmin;};
  double xmax() const {return _xmax;};
  double ymax() const {return _ymax;};
  double zmax() const {return _zmax;};


  // node access
  const Mesh_Node node(unsigned int i) const {return _nodes[i];};

  // property interaction
  void add_phys_property(std::string property_name, double * property_vals);
  void reset_all_properties(std::vector<double> properties);
  void reset_property(std::string property_name, double reset_val);
  //float * get_phys_property_ptr(std::string property_name);


  // grid generation and refinement
  static Static_Mesh * create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static Static_Mesh * create_regular_grid(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);
  //static Mesh * create_unstructured_tri_simple();

  // reading and writing files
  static Static_Mesh * read_MSH(std::string filename, unsigned int byte_offset=0);
  static Static_Mesh * read_NEU(std::string filename, unsigned int byte_offset=0);
  static Static_Mesh * read_CAS(std::string filename, unsigned int byte_offset=0);
  void write_MSH(std::string filename) const;
  void write_NEU(std::string filename) const;
  void write_CAS(std::string filename) const;

private:
  // metadata
  MeshType _mesh_type;
  unsigned int _num_dims, _nodecount, _elementcount;
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  // nodes and elements
  Mesh_Node * _nodes; // array of nodes

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
