/************************************************************************************//**
 * \file function_approx.cpp
 * 
 * File filled with necessary code for function approximation
 *
 ***************************************************************************************/

#include "Mesh.hpp"

//#define _TEST_

using namespace std;

MeshNode::MeshNode(){
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _boundary = false;
  _core_group = 0;
  _num_connections = 0;
}

MeshNode::MeshNode(double x, double y, double z, bool boundary, unsigned int num_connections, unsigned int core_group){
  _x = x;
  _y = y;
  _z = z;
  _boundary = boundary;
  _num_connections = num_connections;
  _core_group = core_group;
}

MeshNode::~MeshNode(){

}

void MeshNode::print_summary() const {
  cout << "  MeshNode: (" << _x << ", " << _y << ", " << _z 
        << ")    boundary?: " << (_boundary? "yes" : "no") 
        << "   num_connections: " << _num_connections
        << "   core_group: " << _core_group << endl;
}

void MeshNode::print_detailed() const{
  cout << "  MeshNode: (" << _x << ", " << _y << ", " << _z 
        << ")    boundary?: " << (_boundary? "yes" : "no") 
        << "   num_connections: " << _num_connections
        << "   core_group: " << _core_group << endl;
}

double MeshNode::dist_sq(const MeshNode & node1, const MeshNode & node2){
  return (node1.x()-node2.x())*(node1.x()-node2.x()) + 
          (node1.y()-node2.y())*(node1.y()-node2.y()) + 
          (node1.z()-node2.z())*(node1.z()-node2.z());
}

double MeshNode::area_tri(const MeshNode & node1, const MeshNode & node2, const MeshNode & node3){
  return fabs((node2.x()*node1.y() - node1.x()*node2.y()) + 
              (node3.x()*node2.y() - node2.x()*node3.y()) + 
              (node1.x()*node3.y() - node3.x()*node1.y()));
}
//******************************************************************************************************

MeshElement::MeshElement(){
  _element_type = EMPTY;
}

MeshElement::MeshElement(std::vector<unsigned int> vertex_inds){
  _vertex_inds = vertex_inds;
  recalc_type();
}

MeshElement::~MeshElement(){

}

void MeshElement::print_summary() const{
  cout << "  MeshElement: num_elements: " << _vertex_inds.size() << " vertices: ";
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    cout << _vertex_inds.at(i) << "  ";
  }
  cout << endl;
  return;
}

void MeshElement::print_detailed() const{
  cout << " Mesh Element: " << _vertex_inds.size() << " vertices: ";
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    cout << _vertex_inds.at(i) << ", ";
  }
  cout << endl;
  return;
}

void MeshElement::remove_vertex(unsigned int vert_ind){
  for (unsigned int i=0; i<_vertex_inds.size(); i++){
    if (_vertex_inds.at(i) == vert_ind){
      // remove index
      _vertex_inds.erase(_vertex_inds.begin()+i);
      break;
    }
  }
  recalc_type();
}

void MeshElement::add_vertex(unsigned int vert_ind, int position){
  if (position == -1){
    _vertex_inds.push_back(vert_ind);
  }
  else{
    _vertex_inds.insert(_vertex_inds.begin() + position, vert_ind);
  }
  recalc_type();
}
 
void MeshElement::recalc_type(){
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
Mesh::Mesh(){
  _xmin = 0;
  _xmax = 0;
  _ymin = 0;
  _ymax = 0;
  _zmin = 0;
  _zmax = 0;
  _mesh_type = MESH_UNKNOWN;
  _num_dims = 0;
}

Mesh::~Mesh(){

}

void Mesh::print_summary() const{
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
  cout << "  Extra Fields (" << _phys_property_names.size() << "): " << endl;
  for (auto i=0; i<_phys_property_names.size(); i++) cout << "      " << _phys_property_names.at(i) << endl;
  cout << "************************************* " << endl;
  cout << " " << endl;

  return;
}

void Mesh::print_detailed() const{

}


const double & Mesh::x() {
  if (_x.size() != _nodes.size()){
    _x.clear();
    _x.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _x.at(i) = _nodes.at(i).x();
    }
  }
  return _x.front();
}


const double & Mesh::y(){
  if (_y.size() != _nodes.size()){
    _y.clear();
    _y.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _y.at(i) = _nodes.at(i).y();
    }
  }
  return _y.front();
}

const double & Mesh::z(){
  if (_z.size() != _nodes.size()){
    _z.clear();
    _z.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _z.at(i) = _nodes.at(i).z();
    }
  }
  return _z.front();
}


const bool & Mesh::boundary(){
  if (_boundary.size() != _nodes.size()){
    _boundary.clear();
    _boundary.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _boundary.at(i) = _nodes.at(i).boundary();
    }
  }
  return _boundary.front();
}


const unsigned int & Mesh::core_group(){
  if (_core_group.size() != _nodes.size()){
    _core_group.clear();
    _core_group.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _core_group.at(i) = _nodes.at(i).core_group();
    }
  }
  return _core_group.front();
}

const unsigned int & Mesh::num_connections(){
  if (_num_connections.size() != _nodes.size()){
    _num_connections.clear();
    _num_connections.resize(_nodes.size());
    for (unsigned int i=0; i<_nodes.size(); i++){
      _num_connections.at(i) = _nodes.at(i).num_connections();
    }
  }
  return _num_connections.front();
}

const double & Mesh::data(std::string fieldname) const{

  if (_phys_properties.at(fieldname).size() != _nodes.size()){
    cout << "physical properties size doesn't match nodes size!" << endl;
    throw -1;
  }

  return _phys_properties.at(fieldname).front();
}

void Mesh::set_nodecount(unsigned int count){
  if (_nodes.size() > 0){
    cout << "cannot expand mesh if it is already allocated!" << endl;
    throw -1;
  }

  _nodes.resize(count);

}

void Mesh::set_elementcount(unsigned int count){
  if (_nodes.size() > 0){
    cout << "cannot expand mesh if it is already allocated!" << endl;
    throw -1;
  }

  _nodes.resize(count);

}

void Mesh::add_phys_property(std::string property_name, const double * property_vals){
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

void Mesh::add_phys_property(std::string property_name, double init_val){
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

void Mesh::reset_property(std::string property_name, double reset_val){
  for (unsigned int i=0; i<_phys_property_names.size(); i++){
    if (_phys_property_names.at(i).compare(property_name) == 0){
      cout << property_name << " is already in the mesh... doing nothing" << endl;
      return;
    }
  }

  return;
}


void Mesh::calc_extents(){
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

Mesh Mesh::read_MSH(string filename, unsigned int byte_offset){
  Mesh outmesh;
  outmesh.read_MSH_internal(filename, byte_offset);

  return outmesh;
}

void Mesh::read_MSH_internal(string filename, unsigned int byte_offset){

}






//*************************************************************************************************************

#ifdef _TEST_

// to compile: g++ -std=c++11 Mesh.cpp -o mesh_test

int main(int argc, char * argv[]){

  // test constructor
  cout << "testing mesh constructor..." << flush;
  Mutable_Mesh * mymesh = new Mutable_Mesh();
  //Mesh * mymesh = new Mesh();
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
  Mesh * mesh_reg_1ds = Mesh::create_regular_grid_n(0.1, 100, 50, 10);
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

