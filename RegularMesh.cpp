#include "RegularMesh.hpp"

using namespace std;

RegularMesh RegularMesh::create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                    unsigned int num_nodes_z){
  RegularMesh mesh_out;
  mesh_out.create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, 0.0, 0.0, 0.0);

  return mesh_out;
}


RegularMesh RegularMesh::create_regular_grid_b(double res, double xmin, double xmax, double ymin, double ymax,
                    double zmin, double zmax){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  RegularMesh mesh_out;

  num_nodes_x = (unsigned int)((xmax-xmin)/res) + 1;
  num_nodes_y = (unsigned int)((ymax-ymin)/res) + 1;
  num_nodes_z = (unsigned int)((zmax-zmin)/res) + 1;

  //cout << "nx: " << num_nodes_x << "  ny: " << num_nodes_y << "  nz: " << num_nodes_z << endl;
  double xcen, ycen, zcen;
  xcen = (xmin + xmax)/2.0;
  ycen = (ymin + ymax)/2.0;
  zcen = (zmin + zmax)/2.0;
  mesh_out.create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, xcen, ycen, zcen);

  return mesh_out;
}



void RegularMesh::create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
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

unsigned int RegularMesh::reg_inds_to_glob_ind(unsigned int i, unsigned int j, unsigned int k){
  return k*(_num_nodes_x*_num_nodes_y) + j*(_num_nodes_x) + i;
}

/*
// regular mesh
class RegularMesh : public Mesh{
public:

  RegularMesh();
  RegularMesh(const RegularMesh & rmesh);
  ~RegularMesh();

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


  // grid generation and refinement
  static RegularMesh * create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static RegularMesh * create_regular_grid_b(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);

private:

  // metadata for regular mesh
  unsigned int _num_nodes_x, _num_nodes_y, _num_nodes_z;
  double _res;
};
*/

