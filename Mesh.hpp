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


enum MeshType{REGULAR=0, UNSTRUCTURED_TRI=1, UNSTRUCTURED_QUAD=2, MESH_UNKNOWN};
enum ElementType{EMPTY, POINT, LINE, TRIANGLE, QUADRANGLE, TETRAHEDRON, HEXAHEDRON, PRISM, PYRAMID, UNKNOWN};


// forward declarations
class MeshNode;
class MeshElement;


class Mesh{
public:
  Mesh();   // constructor
  Mesh(const Mesh & mesh);  // copy constructor
  ~Mesh();    // destructor

  void print_summary() const;
  void print_detailed() const;
  

  // metadata access
  MeshType mesh_type() const {return _mesh_type;};
  unsigned int num_dims() const {return _num_dims;};
  unsigned int nodecount() const {return _nodes.size();};
  unsigned int elementcount() const {return _elements.size();};
  double xmin() const {return _xmin;};
  double ymin() const {return _ymin;};
  double zmin() const {return _zmin;};
  double xmax() const {return _xmax;};
  double ymax() const {return _ymax;};
  double zmax() const {return _zmax;};


  // node and element access
  MeshNode & node(unsigned int i) {return _nodes.at(i);};
  MeshElement & element(unsigned int i) {return _elements.at(i);};
  const MeshNode & node(unsigned int i) const {return _nodes.at(i);};
  const MeshElement & element(unsigned int i) const {return _elements.at(i);};
  
  // property interaction and access
  const double & x();
  const double & y();
  const double & z();
  const bool & boundary();
  const unsigned int & core_group();
  const unsigned int & num_connections();
  const double & data(std::string fieldname) const;

  /*
  // specialty stuff for regular grids
  MeshNode & regular_node(unsigned int i, unsigned int j=0, unsigned int k=0);
  unsigned int reg_num_nodes_x() const {return _num_nodes_x;};
  unsigned int reg_num_nodes_y() const {return _num_nodes_y;};
  unsigned int reg_num_nodes_z() const {return _num_nodes_z;};
  //unsigned int & reg_nodes_boundary(BoundaryLocation loc);
  std::vector<unsigned int> reg_left_inds();
  std::vector<unsigned int> reg_right_inds();
  std::vector<unsigned int> reg_top_inds();
  std::vector<unsigned int> reg_bottom_inds();
  std::vector<unsigned int> reg_front_inds();
  std::vector<unsigned int> reg_back_inds();
  unsigned int reg_inds_to_glob_ind(unsigned int i, unsigned int j=0, unsigned int k=0);
  */

  void add_phys_property(std::string property_name, const double * property_vals);
  void add_phys_property(std::string proprety_name, double init_val);
  void reset_property(std::string property_name, double reset_val=0.0);
  void set_phys_property(std::string property_name, unsigned int i, double val){_phys_properties.at(property_name).at(i) = val;};


  /*
  // grid generation and refinement
  static Mesh * create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static Mesh * create_regular_grid_b(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);
  */

  //static Mesh * create_unstructured_tri_simple();

  // reading and writing files
  static Mesh read_MSH(std::string filename, unsigned int byte_offset=0);
  static Mesh read_XML(std::string filename, unsigned int byte_offset=0);
  static Mesh read_NEU(std::string filename, unsigned int byte_offset=0);
  static Mesh read_CAS(std::string filename, unsigned int byte_offset=0);
  void write_MSH(std::string filename) const;
  void write_NEU(std::string filename) const;
  void write_CAS(std::string filename) const;

//private:
protected:
  // metadata
  MeshType _mesh_type;
  unsigned int _num_dims;
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

  
  // metadata for regular mesh
  unsigned int _num_nodes_x, _num_nodes_y, _num_nodes_z;
  double _res;
  

  // nodes and elements
  std::vector<MeshNode> _nodes; // array of nodes
  std::vector<MeshElement> _elements; // array of elements

  // user-defined propertie for the mesh
  // DYLAN_TODO: RENAME THIS NodeData and ElementData
  std::vector<std::string> _phys_property_names; // the name position in this vector corresponds with the position of the property value in the node
  std::map<std::string, std::vector<double>> _phys_properties;

  // other properties conveniently placed in arrays (on demand) in order to return data
  std::vector<bool> _boundary;
  std::vector<double> _x, _y, _z;
  std::vector<unsigned int> _core_group, _num_connections;

  /*
  void create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen=0.0, double ycen=0.0, double zcen=0.0);
  */

  void calc_extents();

};





// cell class for finite volume methods
class MeshCell{
public:
  unsigned int num_vertices;
  std::vector<unsigned int> vertices;

protected:

private:

};





class MeshElement{
public:
  MeshElement();
  MeshElement(std::vector<unsigned int> vertex_inds);
  ~MeshElement();

  // utils
  void print_summary() const;
  void print_detailed() const;
  double area() const;
  double perimeter() const;

  // member data access
  unsigned int num_vertices() const {return _vertex_inds.size();};
  unsigned int vertex_ind(unsigned int i) const {return _vertex_inds.at(i);};
  //std::vector<unsigned int> vertex_inds() const {return _vertex_inds;};
  const unsigned int & vertex_inds() const {return _vertex_inds.front();};

  // mutators
  void remove_vertex(unsigned int vert_ind);
  void add_vertex(unsigned int vert_ind, int position=-1);
  void set_vertex_inds(std::vector<unsigned int> vertex_inds) {_vertex_inds = vertex_inds;};
  void set_element_type(ElementType type) {_element_type = type;};

private:
  std::vector<unsigned int> _vertex_inds;
  ElementType _element_type; // 

  void recalc_type();

};





class MeshNode{
public:
  MeshNode();
  MeshNode(double x, double y, double z=0.0, bool boundary=false, unsigned int num_connections=0, unsigned int core_group=false);
  ~MeshNode();

  // utils
  void print_summary() const;
  void print_detailed() const;
  static double dist_sq(const MeshNode & node1, const MeshNode & node2);
  static double dist(const MeshNode & node1, const MeshNode & node2) {return sqrt(dist_sq(node1, node2));};
  static double area_tri(const MeshNode & node1, const MeshNode & node2, const MeshNode & node3);

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

// regular mesh

// mutable mesh

// affine cell transformation here


#endif
