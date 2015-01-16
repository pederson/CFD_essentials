/************************************************************************************//**
 * \file function_approx.cpp
 * 
 * File filled with necessary code for function approximation
 *
 ***************************************************************************************/

#include "mesh_class.hpp"

//#define _TEST_

using namespace std;

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

Mesh::Mesh(){
  xmin = 0;
  xmax = 0;
  ymin = 0;
  ymax = 0;
  zmin = 0;
  zmax = 0;
}

Mesh::~Mesh(){
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
void Mesh::calc_extents(){
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

void Mesh::set_mesh_type(MeshType type){
  mesh_type = type;
  return;
}

void Mesh::set_num_dims(unsigned int ndims){
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
void Mesh::set_num_nodes(unsigned int number_of_nodes){
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

void Mesh::set_xmin(double x_min){
  xmin = x_min;
  return;
}

void Mesh::set_ymin(double y_min){
  ymin = y_min;
  return;
}

void Mesh::set_zmin(double z_min){
  zmin = z_min;
  return;
}

void Mesh::set_xmax(double x_max){
  xmax = x_max;
  return;
}

void Mesh::set_ymax(double y_max){
  ymax = y_max;
  return;
}

void Mesh::set_zmax(double z_max){
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
Node * Mesh::get_node_ptr(unsigned int key){
  return mesh_nodes.at(key);
}

/************************************************************************************//**
 * \brief Return the key to a node, given its index
 * 
 *
 *  \param i : index of key within all keys
 *
 ***************************************************************************************/
unsigned int Mesh::get_node_key(unsigned int i){
  return node_keys.at(i);
}

/************************************************************************************//**
 * \brief Insert a new node and automatically generate a new key for it
 * 
 *
 *  \param new_node : node to insert (can be a "new Node()")
 *
 ***************************************************************************************/
void Mesh::add_node(Node * new_node){
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
void Mesh::remove_node(unsigned int key){
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

void Mesh::add_phys_property(string property_name){
  // add property to the list
  phys_property_names.push_back(property_name);

  // make a placeholder in all of the nodes for the new property
  for (auto i=0; i<node_keys.size(); i++){
    mesh_nodes.at(node_keys.at(i))->phys_properties.push_back(0.0);
  }

  return;
}

void Mesh::set_background_properties(vector<double> properties){
  for (auto i=0; i<mesh_nodes.size(); i++){
    mesh_nodes.at(node_keys.at(i))->phys_properties = properties;
  }

  return;
}

unsigned int Mesh::get_phys_property_position(string property_name){
  for (auto i=0; i<phys_property_names.size(); i++){
    if (property_name.compare(phys_property_names.at(i))==0) return i;
  }

  cout << property_name << " is not an available property!" << endl;
  throw -1;
}

float * Mesh::get_phys_property_ptr(string property_name){
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
Mesh * Mesh::create_regular_grid(double res, double xmin, double xmax, double ymin, double ymax,
                      double zmin, double zmax){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  Mesh * mesh_out;

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

int main(int argc, char * argv[]){

  // test constructor
  cout << "testing mesh constructor..." << flush;
  Mesh * mymesh = new Mesh();
  cout << "succeeded" << endl;

  // test num nodes setting
  cout << "testing node setting..." << flush;
  mymesh->set_num_nodes(100);
  cout << "succeeded" << endl;

  // test creation of a regular mesh
  cout << "testing regular grid creation..." << flush;
  Mesh * mesh_reg_1d = Mesh::create_regular_grid(0.1, 100);
  cout << "succeeded" << endl;
  mesh_reg_1d->print_summary();
  mesh_reg_1d->get_node_ptr(50)->print_summary();

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

