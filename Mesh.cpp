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

double Mesh_Node::dist_sq(const Mesh_Node & node1, const Mesh_Node & node2){
  return (node1.x()-node2.x())*(node1.x()-node2.x()) + 
          (node1.y()-node2.y())*(node1.y()-node2.y()) + 
          (node1.z()-node2.z())*(node1.z()-node2.z());
}

double Mesh_Node::area_tri(const Mesh_Node & node1, const Mesh_Node & node2, const Mesh_Node & node3){
  return fabs((node2.x()*node1.y() - node1.x()*node2.y()) + 
              (node3.x()*node2.y() - node2.x()*node3.y()) + 
              (node1.x()*node3.y() - node3.x()*node1.y()));
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
  cout << "  Extra Fields (" << _phys_property_names.size() << "): " << endl;
  for (auto i=0; i<_phys_property_names.size(); i++) cout << "      " << _phys_property_names.at(i) << endl;
  cout << "************************************* " << endl;
  cout << " " << endl;

  return;
}

void Static_Mesh::print_detailed() const{

}

Mesh_Node & Static_Mesh::regular_node(unsigned int i, unsigned int j, unsigned int k){
  return _nodes.at(reg_inds_to_glob_ind(i,j,k));
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
Static_Mesh Static_Mesh::read_MSH(std::string filename, unsigned int byte_offset=0){
  
}

Static_Mesh Static_Mesh::read_NEU(string filename, unsigned int byte_offset=0);
Static_Mesh Static_Mesh::read_CAS(string filename, unsigned int byte_offset=0);
void Static_Mesh::write_MSH(string filename) const;
void Static_Mesh::write_NEU(string filename) const;
void Static_Mesh::write_CAS(string filename) const;
*/

void Static_Mesh::create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen, double ycen, double zcen){
  // declare vars

  // do some input checking

  // initialize things and fill in metadata
  _res = res;
  _num_nodes_x = num_nodes_x;
  _num_nodes_y = num_nodes_y;
  _num_nodes_z = num_nodes_z;
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
        //glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
        glidx = reg_inds_to_glob_ind(k,j,i);
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

unsigned int Static_Mesh::reg_inds_to_glob_ind(unsigned int i, unsigned int j, unsigned int k){
  return k*(_num_nodes_x*_num_nodes_y) + j*(_num_nodes_x) + i;
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

