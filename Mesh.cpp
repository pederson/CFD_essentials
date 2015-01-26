/************************************************************************************//**
 * \file function_approx.cpp
 * 
 * File filled with necessary code for function approximation
 *
 ***************************************************************************************/

#include "Mesh.hpp"

//#define _TEST_

using namespace std;

Mesh_Node::Mesh_Node(){
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _boundary = false;
  _core_group = 0;
  _num_connections = 0;
}

Mesh_Node::Mesh_Node(double x, double y, double z, bool boundary, unsigned int num_connections, unsigned int core_group){
  _x = x;
  _y = y;
  _z = z;
  _boundary = boundary;
  _num_connections = num_connections;
  _core_group = core_group;
}

Mesh_Node::~Mesh_Node(){

}

void Mesh_Node::print_summary() const {
  cout << "  Mesh_Node: (" << _x << ", " << _y << ", " << _z 
        << ")    boundary?: " << (_boundary? "yes" : "no") 
        << "   num_connections: " << _num_connections
        << "   core_group: " << _core_group << endl;
}

void Mesh_Node::print_detailed() const{
  cout << "  Mesh_Node: (" << _x << ", " << _y << ", " << _z 
        << ")    boundary?: " << (_boundary? "yes" : "no") 
        << "   num_connections: " << _num_connections
        << "   core_group: " << _core_group << endl;
}

double Mesh_Node::dist_sq(Mesh_Node * node1, Mesh_Node * node2){
  return (node1->x()-node2->x())*(node1->x()-node2->x()) + 
          (node1->y()-node2->y())*(node1->y()-node2->y()) + 
          (node1->z()-node2->z())*(node1->z()-node2->z());
}

double Mesh_Node::area_tri(Mesh_Node * node1, Mesh_Node * node2, Mesh_Node * node3){
  return fabs((node2->x()*node1->y() - node1->x()*node2->y()) + 
              (node3->x()*node2->y() - node2->x()*node3->y()) + 
              (node1->x()*node3->y() - node3->x()*node1->y()));
}
//******************************************************************************************************

Mesh_Element::Mesh_Element(){
  _element_type = EMPTY;
}

Mesh_Element::Mesh_Element(std::vector<unsigned int> vertex_inds){
  _vertex_inds = vertex_inds;
  recalc_type();
}

Mesh_Element::~Mesh_Element(){

}

void Mesh_Element::print_summary() const{
  cout << "  Mesh_Element: num_elements: " << _vertex_inds.size() << " vertices: ";
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    cout << _vertex_inds.at(i) << "  ";
  }
  cout << endl;
  return;
}

void Mesh_Element::print_detailed() const{
  cout << " Mesh Element: " << _vertex_inds.size() << " vertices: ";
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    cout << _vertex_inds.at(i) << ", ";
  }
  cout << endl;
  return;
}

void Mesh_Element::remove_vertex(unsigned int vert_ind){
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    if (_vertex_inds.at(i) == vert_ind){
      // remove index
      _vertex_inds.erase(_vertex_inds.begin()+i);
      break;
    }
  }
  recalc_type();
}

void Mesh_Element::add_vertex(unsigned int vert_ind, int position){
  if (position == -1){
    _vertex_inds.push_back(vert_ind);
  }
  else{
    _vertex_inds.insert(_vertex_inds.begin() + position, vert_ind);
  }
  recalc_type();
}
 
void Mesh_Element::recalc_type(){
  cout << "doing a recalc doesn't have a defined behavior yet!" << endl;
  switch (_vertex_inds.size()){
    case 0:
      _element_type = EMPTY;
      break;

    case 1:
      _element_type = POINT;
      break;

    case 2:
      _element_type = LINE;
      break;

    case 3:
      _element_type = TRIANGLE;
      break;

    case 4:
      _element_type = QUADRANGLE;
      break;

    case 5:
      _element_type = UNKNOWN;
      break;

    otherwise:
      _element_type = UNKNOWN;
  }
}

//******************************************************************************************************
Static_Mesh::Static_Mesh(){
  _xmin = 0;
  _xmax = 0;
  _ymin = 0;
  _ymax = 0;
  _zmin = 0;
  _zmax = 0;
  _mesh_type = MESH_UNKNOWN;
  _num_dims = 0;
}

Static_Mesh::~Static_Mesh(){

}

void Static_Mesh::print_summary() const{
  cout << " " << endl;
  cout << "******** Static Mesh Summary ******** " << endl;
  if (_nodes.size() == 0){
    cout << "  Mesh is empty!" << endl;
    return;
  }

  cout << "  type: ";
  if (_mesh_type==REGULAR) cout << "REGULAR" << endl;
  else if (_mesh_type==UNSTRUCTURED_TRI) cout << "UNSTRUCTURED TRI" << endl;
  else if (_mesh_type==UNSTRUCTURED_QUAD) cout << "UNSTRUCTURED QUAD" << endl;
  else cout << "MESH_UNKNOWN" << endl;
  cout << "  num_dims: " << _num_dims << endl;
  cout << "  last node: ";
  _nodes.back().print_summary();
  cout << "  last element: ";
  _elements.back().print_summary();
  cout << "  num_nodes: " << _nodes.size() << endl;
  cout << "  num_elements: " << _elements.size() << endl;
  cout << "  x extents: [" << _xmin << ", " << _xmax << "]" << endl;
  cout << "  y extents: [" << _ymin << ", " << _ymax << "]" << endl;
  cout << "  z extents: [" << _zmin << ", " << _zmax << "]" << endl;
  cout << "************************************* " << endl;
  cout << " " << endl;

  return;
}

void Static_Mesh::print_detailed() const{

}

const double & Static_Mesh::x(){
  if (_x.size() != _nodes.size()){
    _x.clear();
    _x.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _x.at(i) = _nodes.at(i).x();
    }
  }
  return _x.front();
}


const double & Static_Mesh::y(){
  if (_y.size() != _nodes.size()){
    _y.clear();
    _y.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _y.at(i) = _nodes.at(i).y();
    }
  }
  return _y.front();
}

const double & Static_Mesh::z(){
  if (_z.size() != _nodes.size()){
    _z.clear();
    _z.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _z.at(i) = _nodes.at(i).z();
    }
  }
  return _z.front();
}


const bool & Static_Mesh::boundary(){
  if (_boundary.size() != _nodes.size()){
    _boundary.clear();
    _boundary.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _boundary.at(i) = _nodes.at(i).boundary();
    }
  }
  return _boundary.front();
}


const unsigned int & Static_Mesh::core_group(){
  if (_core_group.size() != _nodes.size()){
    _core_group.clear();
    _core_group.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _core_group.at(i) = _nodes.at(i).core_group();
    }
  }
  return _core_group.front();
}

const unsigned int & Static_Mesh::num_connections(){
  if (_num_connections.size() != _nodes.size()){
    _num_connections.clear();
    _num_connections.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _num_connections.at(i) = _nodes.at(i).num_connections();
    }
  }
  return _num_connections.front();
}

const double & Static_Mesh::data(std::string fieldname) const{

  if (_phys_properties.at(fieldname).size() != _nodes.size()){
    cout << "physical properties size doesn't match nodes size!" << endl;
    throw -1;
  }

  return _phys_properties.at(fieldname).front();
}

void Static_Mesh::add_phys_property(std::string property_name, const double * property_vals){
  for (unsigned int i=0; i<_phys_property_names.size(); i++){
    if (_phys_property_names.at(i).compare(property_name) == 0){
      cout << property_name << " is already in the mesh... doing nothing" << endl;
      return;
    }
  }

  vector<double> prop;
  //const double * valptr = &property_vals;
  prop.assign(_nodes.size(), 0.0);
  for (unsigned int i=0; i<_nodes.size(); i++) prop.at(i) = property_vals[i];
  _phys_property_names.push_back(property_name);
  _phys_properties[property_name] = prop;
  return;
}

void Static_Mesh::add_phys_property(std::string property_name, double init_val){
  for (unsigned int i=0; i<_phys_property_names.size(); i++){
    if (_phys_property_names.at(i).compare(property_name) == 0){
      cout << property_name << " is already in the mesh... doing nothing" << endl;
      return;
    }
  }

  vector<double> prop;
  prop.assign(_nodes.size(), init_val);
  _phys_property_names.push_back(property_name);
  _phys_properties[property_name] = prop;
  return;
}

void Static_Mesh::reset_property(std::string property_name, double reset_val){
  for (unsigned int i=0; i<_phys_property_names.size(); i++){
    if (_phys_property_names.at(i).compare(property_name) == 0){
      cout << property_name << " is already in the mesh... doing nothing" << endl;
      return;
    }
  }

  return;
}


Static_Mesh * Static_Mesh::create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                    unsigned int num_nodes_z){
  Static_Mesh * mesh_out = new Static_Mesh();
  mesh_out->create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, 0.0, 0.0, 0.0);

  return mesh_out;
}


Static_Mesh * Static_Mesh::create_regular_grid_b(double res, double xmin, double xmax, double ymin, double ymax,
                    double zmin, double zmax){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  Static_Mesh * mesh_out = new Static_Mesh();

  num_nodes_x = (unsigned int)((xmax-xmin)/res) + 1;
  num_nodes_y = (unsigned int)((ymax-ymin)/res) + 1;
  num_nodes_z = (unsigned int)((zmax-zmin)/res) + 1;

  //cout << "nx: " << num_nodes_x << "  ny: " << num_nodes_y << "  nz: " << num_nodes_z << endl;
  double xcen, ycen, zcen;
  xcen = (xmin + xmax)/2.0;
  ycen = (ymin + ymax)/2.0;
  zcen = (zmin + zmax)/2.0;
  mesh_out->create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, xcen, ycen, zcen);

  return mesh_out;
}


/*
Static_Mesh * Static_Mesh::read_MSH(std::string filename, unsigned int byte_offset=0){
  
}

Static_Mesh * Static_Mesh::read_NEU(std::string filename, unsigned int byte_offset=0);
Static_Mesh * Static_Mesh::read_CAS(std::string filename, unsigned int byte_offset=0);
void Static_Mesh::write_MSH(std::string filename) const;
void Static_Mesh::write_NEU(std::string filename) const;
void Static_Mesh::write_CAS(std::string filename) const;
*/

void Static_Mesh::create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen, double ycen, double zcen){
  // declare vars

  // do some input checking

  // initialize things and fill in metadata
  _mesh_type = REGULAR;
  _nodes.resize(num_nodes_x*num_nodes_y*num_nodes_z);
  if (num_nodes_y==1 && num_nodes_z==1){
    _num_dims = 1;
  }
  else if (num_nodes_z==1){
    _num_dims = 2;
  }
  else{
    _num_dims = 3;
  }

  // create nodes
  unsigned int glidx, i, j, k;
  for (i=0; i<num_nodes_z; i++){
    for (j=0; j<num_nodes_y; j++){
      for (k=0; k<num_nodes_x; k++){
        glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
        _nodes.at(glidx).set_z(double(i)*res);
        _nodes.at(glidx).set_y(double(j)*res);
        _nodes.at(glidx).set_x(double(k)*res);
      }
    }
  }
  // set boundary properties for nodes
  if (num_nodes_y==1 && num_nodes_z==1){
    i=0; j=0;
    for (k=0; k<num_nodes_x; k++){
      _nodes.at(k).set_num_connections(2);
    }
    _nodes.at(0).set_boundary(true);
    _nodes.at(0).set_num_connections(1);
    _nodes.at(num_nodes_x-1).set_boundary(true);
    _nodes.at(num_nodes_x-1).set_num_connections(1);
  }
  else if (num_nodes_z==1){
    i=0;
    for (j=0; j<num_nodes_y; j++){
      for (k=0; k<num_nodes_x; k++){
        glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
        if (k==0 || k==num_nodes_x-1 || j==0 || j==num_nodes_y-1){
          _nodes.at(glidx).set_boundary(true);
          _nodes.at(glidx).set_num_connections(3);
        }
      }
    }
    _nodes.at(0).set_num_connections(2);
    _nodes.at(num_nodes_x-1).set_num_connections(2);
    _nodes.at((num_nodes_y-1)*num_nodes_x).set_num_connections(2);
    _nodes.at((num_nodes_y-1)*num_nodes_x+num_nodes_x-1).set_num_connections(2);
  }
  else{
    for (i=0; i<num_nodes_z; i++){
      for (j=0; j<num_nodes_y; j++){
        for (k=0; k<num_nodes_x; k++){
          glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
          if (k==0 || k==num_nodes_x-1 || j==0 || j==num_nodes_y-1 || i==0 || i==num_nodes_z-1){
            _nodes.at(glidx).set_boundary(true);
            _nodes.at(glidx).set_num_connections(5);
          }
        }
      }
    }
    _nodes.at(0).set_num_connections(3);
    _nodes.at(num_nodes_x-1).set_num_connections(3);
    _nodes.at((num_nodes_y)*(num_nodes_x)-1).set_num_connections(3);
    _nodes.at((num_nodes_y-1)*(num_nodes_x)+1).set_num_connections(3);

    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1).set_num_connections(3);
    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1 + num_nodes_x-1).set_num_connections(3);
    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1 + (num_nodes_y)*(num_nodes_x)-1).set_num_connections(3);
    _nodes.at((num_nodes_z)*(num_nodes_x)*(num_nodes_y)-1).set_num_connections(3);
  }
  

  // create elements
  unsigned int blf, tlf, brf, trf, blb, tlb, brb, trb, nex, ney, nez; 
  
  if (num_nodes_y==1 && num_nodes_z==1){
    nex = num_nodes_x-1;
    _elements.resize((num_nodes_x-1));
    i=0; j=0;
    for (k=0; k<_elements.size(); k++){
      glidx = i*(nex*ney) + j*(nex) + k;
      blf = k;
      brf = k+1;
      _elements.at(glidx).set_vertex_inds({blf, brf});
      _elements.at(glidx).set_element_type(LINE);
    }
  }
  else if (num_nodes_z==1){
    nex = num_nodes_x-1;
    ney = num_nodes_y-1;
    _elements.resize((num_nodes_x-1)*(num_nodes_y-1));
    for (j=0; j<ney; j++){
      for (k=0; k<nex; k++){
        glidx = i*(nex*ney) + j*(nex) + k;
        blf = (nex+1)*(j) + k;
        brf = (nex+1)*(j) + k+1;
        trf = (nex+1)*(j+1) + k+1;
        tlf = (nex+1)*(j+1) + k;
        _elements.at(glidx).set_vertex_inds({blf, brf, trf, tlf});
        _elements.at(glidx).set_element_type(QUADRANGLE);
      }
    }
  }
  else{
    cout << "made it here" << endl;
    _elements.resize((num_nodes_x-1)*(num_nodes_y-1)*(num_nodes_z-1));
    nex = num_nodes_x-1;
    ney = num_nodes_y-1;
    nez = num_nodes_z-1;
    for (i=0; i<nez; i++){
      for (j=0; j<ney; j++){
        for (k=0; k<nex; k++){
          glidx = i*(nex*ney) + j*(nex) + k;
          blf = (nex+1)*(ney+1)*(i) + (nex+1)*(j) + k;
          brf = (nex+1)*(ney+1)*(i) + (nex+1)*(j) + k+1;
          trf = (nex+1)*(ney+1)*(i) + (nex+1)*(j+1) + k+1;
          tlf = (nex+1)*(ney+1)*(i) + (nex+1)*(j+1) + k;
          blb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j) + k;
          brb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j) + k+1;
          trb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j+1) + k+1;
          tlb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j+1) + k;
          _elements.at(glidx).set_vertex_inds({blf, brf, trf, tlf, blb, brb, trb, tlb});
          _elements.at(glidx).set_element_type(HEXAHEDRON);
        }
      }
    }
  }

  // translate the center to 0,0
  double xc, yc, zc;
  xc = ((num_nodes_x-1)*res)/2.0;
  yc = ((num_nodes_y-1)*res)/2.0;
  zc = ((num_nodes_z-1)*res)/2.0;

  // translate to the correct centerpoint
  for (auto i=0; i<_nodes.size(); i++){
    //nd = _nodes.at(i);
    _nodes.at(i).set_x(_nodes.at(i).x() + xcen - xc);
    _nodes.at(i).set_y(_nodes.at(i).y() + ycen - yc);
    _nodes.at(i).set_z(_nodes.at(i).z() + zcen - zc);
  }
  calc_extents();

  return;
}

void Static_Mesh::calc_extents(){
  // declare vars

  _xmax = _nodes.at(0).x(); _xmin = _nodes.at(0).x(); _ymax = _nodes.at(0).y(); _ymin = _nodes.at(0).y(); _zmax = _nodes.at(0).z(); _zmin = _nodes.at(0).z();
  for (unsigned int i=1; i<_nodes.size(); i++){
    if (_nodes.at(i).x() > _xmax) _xmax = _nodes.at(i).x();
    if (_nodes.at(i).x() < _xmin) _xmin = _nodes.at(i).x();
    if (_nodes.at(i).y() > _ymax) _ymax = _nodes.at(i).y();
    if (_nodes.at(i).y() < _ymin) _ymin = _nodes.at(i).y();
    if (_nodes.at(i).z() > _zmax) _zmax = _nodes.at(i).z();
    if (_nodes.at(i).z() < _zmin) _zmin = _nodes.at(i).z();
  }
  return;
}













//******************************************************************************************************

Node::Node(){
  x = 0;
  y = 0; 
  z = 0;

  boundary = false;
  core_group = 0;
}

Node::~Node(){
  
}

/************************************************************************************//**
 * \brief Print a brief summary of Node values
 * 
 ***************************************************************************************/
void Node::print_summary(){
  cout << "Node Summary: " << endl;
  cout << "  key: " << key << endl;
  cout << "  x: " << x << " y: " << y << " z: " << z << endl;
  cout << "  On boundary?: " << (boundary? "true":"false") << endl;
  cout << "  # of neighbors: " << neighbor_keys.size() << endl;
  for (unsigned int i=0; i<neighbor_keys.size(); i++) cout << "    neighbor key: " << neighbor_keys.at(i) << endl;

  return;
}

//******************************************************************************************************

Mutable_Mesh::Mutable_Mesh(){
  xmin = 0;
  xmax = 0;
  ymin = 0;
  ymax = 0;
  zmin = 0;
  zmax = 0;
}

Mutable_Mesh::~Mutable_Mesh(){
  if (node_keys.size() > 0){
    unsigned int nnodes = node_keys.size();
    for (unsigned int i=0; i<nnodes; i++) {
      delete mesh_nodes.at(node_keys.at(i));
    }
  }
}

/************************************************************************************//**
 * \brief Print a brief summary of Mesh values
 * 
 ***************************************************************************************/
void Mutable_Mesh::print_summary(){
  cout << "Mesh Summary: " << endl;
  if (node_keys.size() == 0){
    cout << "  Mesh is empty!" << endl;
    return;
  }

  cout << "  type: ";
  if (mesh_type==REGULAR) cout << "REGULAR" << endl;
  else if (mesh_type==UNSTRUCTURED_TRI) cout << "UNSTRUCTURED TRI" << endl;
  else if (mesh_type==UNSTRUCTURED_QUAD) cout << "UNSTRUCTURED QUAD" << endl;
  else cout << "UNKNOWN" << endl;
  cout << "  num_dims: " << num_dims << endl;
  cout << "  last key: " << node_keys.back() << endl;
  cout << "  num_nodes: " << node_keys.size() << endl;
  cout << "  x extents: [" << xmin << ", " << xmax << "]" << endl;
  cout << "  y extents: [" << ymin << ", " << ymax << "]" << endl;
  cout << "  z extents: [" << zmin << ", " << zmax << "]" << endl;

  return;
}

/************************************************************************************//**
 * \brief Recalculate x,y,z extents of the node data
 *
 ***************************************************************************************/
void Mutable_Mesh::calc_extents(){
  // declare vars
  Node * this_node;
  this_node = mesh_nodes.at(node_keys.at(0));

  xmax = this_node->x; xmin = this_node->x; ymax = this_node->y; ymin = this_node->y; zmax = this_node->z; zmin = this_node->z;
  for (unsigned int i=1; i<node_keys.size(); i++){
    this_node = mesh_nodes.at(node_keys.at(i));
    if (this_node->x > xmax) xmax = this_node->x;
    if (this_node->x < xmin) xmin = this_node->x;
    if (this_node->y > ymax) ymax = this_node->y;
    if (this_node->y < ymin) ymin = this_node->y;
    if (this_node->z > zmax) zmax = this_node->z;
    if (this_node->z < zmin) zmin = this_node->z;
  }

  return;
}

void Mutable_Mesh::set_mesh_type(MeshType type){
  mesh_type = type;
  return;
}

void Mutable_Mesh::set_num_dims(unsigned int ndims){
  num_dims = ndims;
  return;
}


/************************************************************************************//**
 * \brief Set the number of nodes in the mesh
 * 
 *  If the number of nodes is larger than the current number, it will
 *  allocate new nodes up to the new number. Otherwise it will throw an error
 *
 *  \param number_of_nodes : number of nodes to extend to
 *
 ***************************************************************************************/
void Mutable_Mesh::set_num_nodes(unsigned int number_of_nodes){
  unsigned int new_key, begin_size;

  if (number_of_nodes < node_keys.size()){
    cout << "ERROR: cannot set a number of nodes smaller than current size" << endl;
    throw -1;
  }

  begin_size = node_keys.size();
  
  // generate new keys
  if (begin_size == 0) new_key = 0;
  else new_key = node_keys.back() + 1;
  node_keys.resize(number_of_nodes);

  // create new Node objects and put them in the mesh_nodes list
  for (unsigned int i=begin_size; i<number_of_nodes; i++){
    node_keys.at(i) = new_key;
    mesh_nodes[new_key] = new Node();
    mesh_nodes.at(new_key)->key = new_key;

    new_key++;
  }
  return;
}

/*
void Mesh::set_offsets(double x_off, double y_off, double z_off){
  x_offset = x_off;
  y_offset = y_off;
  z_offset = z_off;
  return;
}
*/

void Mutable_Mesh::set_xmin(double x_min){
  xmin = x_min;
  return;
}

void Mutable_Mesh::set_ymin(double y_min){
  ymin = y_min;
  return;
}

void Mutable_Mesh::set_zmin(double z_min){
  zmin = z_min;
  return;
}

void Mutable_Mesh::set_xmax(double x_max){
  xmax = x_max;
  return;
}

void Mutable_Mesh::set_ymax(double y_max){
  ymax = y_max;
  return;
}

void Mutable_Mesh::set_zmax(double z_max){
  zmax = z_max;
  return;
}

/************************************************************************************//**
 * \brief Return a pointer to node, given a key
 * 
 *
 *  \param key : key to the desired node pointer
 *
 ***************************************************************************************/
Node * Mutable_Mesh::get_node_ptr(unsigned int key){
  return mesh_nodes.at(key);
}

/************************************************************************************//**
 * \brief Return the key to a node, given its index
 * 
 *
 *  \param i : index of key within all keys
 *
 ***************************************************************************************/
unsigned int Mutable_Mesh::get_node_key(unsigned int i){
  return node_keys.at(i);
}

/************************************************************************************//**
 * \brief Insert a new node and automatically generate a new key for it
 * 
 *
 *  \param new_node : node to insert (can be a "new Node()")
 *
 ***************************************************************************************/
void Mutable_Mesh::add_node(Node * new_node){
  unsigned int new_key;
  new_key = node_keys.back()+1;
  node_keys.push_back(new_key);
  mesh_nodes[new_key] = new_node;
  mesh_nodes.at(new_key)->key = new_key;

  // check that the new node doesn't exceed the mins and maxs
  if (new_node->x < xmin) xmin = new_node->x;
  else if (new_node->x > xmax) xmax = new_node->x;
  if (new_node->y < ymin) ymin = new_node->y;
  else if (new_node->y > ymax) ymax = new_node->y;
  if (new_node->z < zmin) zmin = new_node->z;
  else if (new_node->z > zmax) zmax = new_node->z;

  return;
}

/************************************************************************************//**
 * \brief Remove a node and it's neighbor connections if you know its key
 * 
 *
 *  \param key : key of Node to remove
 *
 ***************************************************************************************/
void Mutable_Mesh::remove_node(unsigned int key){
  // declare vars
  unsigned int nneighbs, nneighb_j, i;

  i = key;

  // delete neighbor associations to this point
  nneighbs = mesh_nodes.at(i)->neighbor_keys.size();
  for (unsigned int j=0; j<nneighbs; j++){ // for each neighbor
    // search for this association
    nneighb_j = mesh_nodes.at(mesh_nodes.at(i)->neighbor_keys.at(j))->neighbor_keys.size();
    for (unsigned int k=0; k<nneighb_j; k++){ // search through neighbor keys
      if (mesh_nodes.at(mesh_nodes.at(i)->neighbor_keys.at(j))->neighbor_keys.at(k) == i){
        mesh_nodes.at(mesh_nodes.at(i)->neighbor_keys.at(j))->neighbor_keys.erase(mesh_nodes.at(mesh_nodes.at(i)->neighbor_keys.at(j))->neighbor_keys.begin() + k);
        break;
      }
    }
  }

  // delete the node
  delete mesh_nodes.at(i); // this is commented out temporarily because it is deleted upon destroying the object
  mesh_nodes.erase(i);

  // delete the key
  for (unsigned int j=0; j<node_keys.size(); j++){
    if (node_keys.at(j) == key){
      node_keys.erase(node_keys.begin() + j);
      break;
    }
  }


  return;
}

void Mutable_Mesh::add_phys_property(string property_name){
  // add property to the list
  phys_property_names.push_back(property_name);

  // make a placeholder in all of the nodes for the new property
  for (auto i=0; i<node_keys.size(); i++){
    mesh_nodes.at(node_keys.at(i))->phys_properties.push_back(0.0);
  }

  return;
}

void Mutable_Mesh::set_background_properties(vector<double> properties){
  for (auto i=0; i<mesh_nodes.size(); i++){
    mesh_nodes.at(node_keys.at(i))->phys_properties = properties;
  }

  return;
}

unsigned int Mutable_Mesh::get_phys_property_position(string property_name){
  for (auto i=0; i<phys_property_names.size(); i++){
    if (property_name.compare(phys_property_names.at(i))==0) return i;
  }

  cout << property_name << " is not an available property!" << endl;
  throw -1;
}

float * Mutable_Mesh::get_phys_property_ptr(string property_name){
  float * prop;

  unsigned int pos = get_phys_property_position(property_name);

  prop = new float[mesh_nodes.size()];
  for (auto i=0; i<mesh_nodes.size(); i++){
    prop[i] = mesh_nodes.at(node_keys.at(i))->phys_properties.at(pos);
  }

  return prop;
}


/************************************************************************************//**
 * \brief Create a 1d, 2d, or 3d regular grid
 * 
 *
 *  \param res : resolution of the grid
 *  \param num_nodes_x : number of nodes in x direction
 *  \param num_nodes_y : number of nodes in y direction
 *  \param num_nodes_z : number of nodes in z direction
 *
 ***************************************************************************************/
Mutable_Mesh * Mutable_Mesh::create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z){
  // declare vars
  Mutable_Mesh * mesh_out;
  Node * node_spawn;
  unsigned int nodes_total;
  unsigned int glob_idx;

  // do some input checking

  // fill in the extra member data
  mesh_out = new Mutable_Mesh();
  mesh_out->set_mesh_type(REGULAR);
  nodes_total = num_nodes_x*num_nodes_y*num_nodes_z;
  mesh_out->set_num_nodes(nodes_total);

  // begin creating nodes with x, y, z points
  if (num_nodes_y == 1 && num_nodes_z==1){ // x points only

    mesh_out->set_num_dims(1);
    mesh_out->set_xmax(double(num_nodes_x-1)*res);

    for (unsigned int i=0; i<nodes_total; i++){

      // overall index
      glob_idx = i;
      node_spawn = mesh_out->get_node_ptr(glob_idx);
      node_spawn->x = double(i)*res;

      // set neighbors and boundaries
      if (i == 0 || i == num_nodes_x-1){
        node_spawn->boundary = true;
        if (i==0) node_spawn->neighbor_keys.push_back(i+1);
        else node_spawn->neighbor_keys.push_back(i-1);
      } 
      else{
        node_spawn->boundary = false;
        node_spawn->neighbor_keys.push_back(i+1);
        node_spawn->neighbor_keys.push_back(i-1);
      }
    }

  }
  else if (num_nodes_z == 1){ // x and y points only

    mesh_out->set_num_dims(2);
    mesh_out->set_xmax(double(num_nodes_x-1)*res);
    mesh_out->set_ymax(double(num_nodes_y-1)*res);

    for (unsigned int i=0; i<num_nodes_x; i++){
      for (unsigned int j=0; j<num_nodes_y; j++){

        glob_idx = i*(num_nodes_y) + j;
        node_spawn = mesh_out->get_node_ptr(glob_idx);
        node_spawn->x = double(i)*res;
        node_spawn->y = double(j)*res;

        // boundaries
        if (i == 0 || j == 0 || i == num_nodes_x-1 || j == num_nodes_y-1){
          node_spawn->boundary = true;

          if (i==0){
            if (j==0){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
              node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
            }
            else if(j==num_nodes_y-1){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
              node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
            }
            else{
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
              node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
            }

          }
          else if (i==num_nodes_x-1){
            if (j==0){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
              node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
            }
            else if(j==num_nodes_y-1){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
              node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
            }
            else{
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
              node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
            }

          }
          else{
            if (j==0){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
              node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
              node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
            }
            else if(j==num_nodes_y-1){
              node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
              node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
              node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
            }
          }
        }
        else {
          node_spawn->boundary = false;

          node_spawn->neighbor_keys.push_back(i*num_nodes_y + j+1); // top
          node_spawn->neighbor_keys.push_back(i*num_nodes_y + j-1); // bottom
          node_spawn->neighbor_keys.push_back((i-1)*num_nodes_y + j); // left
          node_spawn->neighbor_keys.push_back((i+1)*num_nodes_y + j); // right
        }
      }
    }

  }
  else{ // x, y, and z points

    mesh_out->set_num_dims(3);
    mesh_out->set_xmax(double(num_nodes_x-1)*res);
    mesh_out->set_ymax(double(num_nodes_y-1)*res);
    mesh_out->set_zmax(double(num_nodes_z-1)*res);

    for (unsigned int i=0; i<num_nodes_x; i++){
      for (unsigned int j=0; j<num_nodes_y; j++){
        for (unsigned int k=0; k<num_nodes_z; k++){

          glob_idx = i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k;
          node_spawn = mesh_out->get_node_ptr(glob_idx);
          node_spawn->x = double(i)*res;
          node_spawn->y = double(j)*res;
          node_spawn->z = double(k)*res;


          // boundaries
        if (i == 0 || j == 0 || k == 0 || i == num_nodes_x-1 || j == num_nodes_y-1 || k == num_nodes_z-1){
          node_spawn->boundary = true;

          if (i==0){
            if (j==0){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }
            else if (j==num_nodes_y-1){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange

              }
            }
            else{
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }
          }

          else if (i==num_nodes_x-1){
            if (j==0){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }
            else if (j==num_nodes_y-1){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange

              }
            }
            else{
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }
          }

          else{
            if (j==0){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }
            else if (j==num_nodes_y-1){
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange

              }
            }
            else{
              if (k==0){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
              }
              else if(k==num_nodes_z-1){
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
              else{
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
                node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
                node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
                node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
              }
            }

          }
        }

        else {
          node_spawn->boundary = false;

          node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j+1)*num_nodes_z + k); // top
          node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + (j-1)*num_nodes_z + k); // bottom
          node_spawn->neighbor_keys.push_back((i-1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // left
          node_spawn->neighbor_keys.push_back((i+1)*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k); // right
          node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k+1); // charm
          node_spawn->neighbor_keys.push_back(i*(num_nodes_y*num_nodes_z) + j*num_nodes_z + k-1); // strange
        }
      }
      }
    }
  }

  return mesh_out;
}

/************************************************************************************//**
 * \brief Create a 1d, 2d, or 3d regular grid, given the bounds
 * 
 *
 *  \param res : resolution of the grid
 *  \param xmin : min x value
 *  \param ymin : min y value
 *  \param zmin : min z value
 *
 ***************************************************************************************/
Mutable_Mesh * Mutable_Mesh::create_regular_grid(double res, double xmin, double xmax, double ymin, double ymax,
                      double zmin, double zmax){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  Mutable_Mesh * mesh_out;

  num_nodes_x = (unsigned int)((xmax-xmin)/res) + 1;
  num_nodes_y = (unsigned int)((ymax-ymin)/res) + 1;
  num_nodes_z = (unsigned int)((zmax-zmin)/res) + 1;

  //cout << "nx: " << num_nodes_x << "  ny: " << num_nodes_y << "  nz: " << num_nodes_z << endl;

  mesh_out = create_regular_grid(res, num_nodes_x, num_nodes_y, num_nodes_z);

  // translate to the correct centerpoint
  Node * nd;
  for (auto i=0; i<mesh_out->get_num_nodes(); i++){
    nd = mesh_out->get_node_ptr(mesh_out->get_node_key(i));

    nd->x += xmin;
    nd->y += ymin;
    nd->z += zmin;
  }
  mesh_out->calc_extents();

  return mesh_out;
}

//*************************************************************************************************************

#ifdef _TEST_

// to compile: g++ -std=c++11 Mesh.cpp -o mesh_test

int main(int argc, char * argv[]){

  // test constructor
  cout << "testing mesh constructor..." << flush;
  Mutable_Mesh * mymesh = new Mutable_Mesh();
  //Static_Mesh * mymesh = new Static_Mesh();
  cout << "succeeded" << endl;

  // test num nodes setting
  cout << "testing node setting..." << flush;
  mymesh->set_num_nodes(100);
  cout << "succeeded" << endl;

  // test creation of a regular mesh
  cout << "testing regular grid creation..." << flush;
  Mutable_Mesh * mesh_reg_1d = Mutable_Mesh::create_regular_grid(0.1, 100);
  cout << "succeeded" << endl;
  mesh_reg_1d->print_summary();
  mesh_reg_1d->get_node_ptr(50)->print_summary();

  // test creation of a regular mesh
  cout << "testing regular grid creation..." << flush;
  Static_Mesh * mesh_reg_1ds = Static_Mesh::create_regular_grid_n(0.1, 100, 50, 10);
  cout << "succeeded" << endl;
  mesh_reg_1ds->print_summary();
  mesh_reg_1ds->node(50).print_summary();
  mesh_reg_1ds->element(50).print_summary();

  // test the node addition
  cout << "testing node insertion..." << flush;
  mesh_reg_1d->add_node(new Node);
  cout << "succeeded" << endl;
  mesh_reg_1d->print_summary();

  // test the node deletion
  cout << "testing node deletion..." << flush;
  mesh_reg_1d->remove_node(100);
  cout << "succeeded" << endl;
  mesh_reg_1d->print_summary();


  cout << "testing mesh deletion..." << flush;
  delete mesh_reg_1d;
  cout << "succeeded" << endl;

  return 0;
}

#endif

