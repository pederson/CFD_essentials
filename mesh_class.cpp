#include "mesh_class.hpp"

using namespace std;

#define _TEST_

//DYLAN_TODO: Use a std::map instead of vectors to associate and add/delete data

Node::Node(){
  x = 0;
  y = 0; 
  z = 0;

}

Node::~Node(){
  
}

//******************************************************************************************************

Mesh::Mesh(){
  x_offset = 0;
  y_offset = 0;
  z_offset = 0;

  xmin = 0;
  xmax = 0;
  ymin = 0;
  ymax = 0;
  zmin = 0;
  zmax = 0;
}

Mesh::~Mesh(){
  if (mesh_nodes.size() > 0){
    unsigned int nnodes = node_keys.size();
    for (unsigned int i=0; i<nnodes; i++) {
      delete mesh_nodes.at(node_keys[i]);
    }
  }
}

void Mesh::print_summary(){
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
  cout << "  num_nodes: " << mesh_nodes.size() << endl;
  cout << "  x extents: [" << xmin + x_offset << ", " << xmax + x_offset << "]" << endl;
  cout << "  y extents: [" << ymin + y_offset << ", " << ymax + y_offset << "]" << endl;
  cout << "  z extents: [" << zmin + z_offset << ", " << zmax + z_offset << "]" << endl;

  return;
}

MeshType Mesh::get_mesh_type(){
  return mesh_type;
}

void Mesh::set_mesh_type(MeshType type){
  mesh_type = type;
  return;
}

unsigned int Mesh::get_num_dims(){
  return num_dims;
}

void Mesh::set_num_dims(unsigned int ndims){
  num_dims = ndims;
  return;
}

unsigned int Mesh::get_num_nodes(){
  return mesh_nodes.size();
}

void Mesh::set_num_nodes(unsigned int number_of_nodes){
  unsigned int begin_size = mesh_nodes.size();
  mesh_nodes.resize(number_of_nodes);

  // create new Node objects and put them in the mesh_nodes list
  for (unsigned int i=begin_size; i<number_of_nodes; i++) mesh_nodes[i] = new Node();
  return;
}

void Mesh::set_offsets(double x_off, double y_off, double z_off){
  x_offset = x_off;
  y_offset = y_off;
  z_offset = z_off;
  return;
}

double Mesh::get_offset_x(){
  return x_offset;
}

double Mesh::get_offset_y(){
  return y_offset;
}

double Mesh::get_offset_z(){
  return z_offset;
}

double Mesh::get_xmin(){
  return xmin;
}

void Mesh::set_xmin(double x_min){
  xmin = x_min;
  return;
}

double Mesh::get_ymin(){
  return ymin;
}

void Mesh::set_ymin(double y_min){
  ymin = y_min;
  return;
}

double Mesh::get_zmin(){
  return zmin;
}

void Mesh::set_zmin(double z_min){
  zmin = z_min;
  return;
}

double Mesh::get_xmax(){
  return xmax;
}

void Mesh::set_xmax(double x_max){
  xmax = x_max;
  return;
}

double Mesh::get_ymax(){
  return ymax;
}

void Mesh::set_ymax(double y_max){
  ymax = y_max;
  return;
}

double Mesh::get_zmax(){
  return zmax;
}

void Mesh::set_zmax(double z_max){
  zmax = z_max;
  return;
}

Node * Mesh::get_node_ptr(unsigned int i){
  return mesh_nodes[i];
}

void Mesh::add_node(Node * new_node){
  mesh_nodes.push_back(new_node);

  // check that the new node doesn't exceed the mins and maxs
  if (new_node->x < xmin) xmin = new_node->x;
  else if (new_node->x > xmax) xmax = new_node->x;
  if (new_node->y < ymin) ymin = new_node->y;
  else if (new_node->y > ymax) ymax = new_node->y;
  if (new_node->z < zmin) zmin = new_node->z;
  else if (new_node->z > zmax) zmax = new_node->z;

  return;
}

//DYLAN_TODO : this doesn't actually delete the nodes... it should
void Mesh::remove_node(unsigned int i){
  // declare vars
  unsigned int nneighbs, nneighb_j;

  // delete neighbor associations to this point
  nneighbs = mesh_nodes[i]->neighbor_index.size();
  for (unsigned int j=0; j<nneighbs; j++){
    // search for this association
    nneighb_j = mesh_nodes[mesh_nodes[i]->neighbor_index[j]]->neighbor_index.size();
    for (unsigned int k=0; k<nneighb_j; k++){
      if (mesh_nodes[mesh_nodes[i]->neighbor_index[j]]->neighbor_index[k] == i){
        mesh_nodes[mesh_nodes[i]->neighbor_index[j]]->neighbor_index.erase(mesh_nodes[mesh_nodes[i]->neighbor_index[j]]->neighbor_index.begin());
        break;
      }
    }
  }

  // delete the node
  //delete mesh_nodes[i]; // this is commented out temporarily because it is deleted upon destroying the object

  return;
}

Mesh * Mesh::create_regular_grid(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z){
  // declare vars
  Mesh * mesh_out;
  Node * node_spawn;
  unsigned int nodes_total;
  unsigned int glob_idx;

  // do some input checking

  // fill in the extra member data
  mesh_out = new Mesh();
  mesh_out->set_mesh_type(REGULAR);
  nodes_total = num_nodes_x*num_nodes_y*num_nodes_z;
  mesh_out->set_num_nodes(nodes_total);

  // begin creating nodes with x, y, z points
  if (num_nodes_y == 1 && num_nodes_z==1){ // x points only

    mesh_out->set_num_dims(1);
    mesh_out->set_xmax(double(num_nodes_x-1)*res);

    for (unsigned int i=0; i<nodes_total; i++){
      glob_idx = i;
      node_spawn = mesh_out->get_node_ptr(glob_idx);
      node_spawn->x = double(i)*res;
      node_spawn->boundary = false;
      
    }

    // boundaries
    node_spawn = mesh_out->get_node_ptr(0);
    node_spawn->boundary = true;
    node_spawn = mesh_out->get_node_ptr(mesh_out->get_num_nodes()-1);
    node_spawn->boundary = true;
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
        if (i == 0 || j == 0 || i == num_nodes_x-1 || j == num_nodes_y-1) node_spawn->boundary = true;
        else node_spawn->boundary = false;
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
          if (i == 0 || j == 0 || k == 0 || i == num_nodes_x-1 || j == num_nodes_y-1 || k == num_nodes_z-1) node_spawn->boundary = true;
          else node_spawn->boundary = false;
        }
      }
    }
  }

  return mesh_out;
}

//*************************************************************************************************************

#ifdef _TEST_

int main(int argc, char * argv[]){

  // test creation of a regular mesh
  Mesh * mesh_reg_1d = Mesh::create_regular_grid(0.1, 100, 50);
  mesh_reg_1d->print_summary();

  // test the node addition
  mesh_reg_1d->add_node(new Node);
  mesh_reg_1d->print_summary();

  // test the node deletion
  mesh_reg_1d->remove_node(100);
  mesh_reg_1d->print_summary();


  cout << "about to delete mesh" << endl;
  delete mesh_reg_1d;

  return 0;
}

#endif

