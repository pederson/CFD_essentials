#include "VisualixerMesh.hpp"

using namespace std;

//#define _TEST_
vector<visualixer*> _vMeshInstances;

//*************************************** Mesh Visualixer Class *******************************************
mesh_visualixer::mesh_visualixer(){

	_vMeshInstances.push_back(this);

	visualixer_active = false;
	window_name = "Mesh Visualixer";
	rotation_lock = false;
	_colorby = nullptr;
	_color_alpha = nullptr;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;

	_num_point_elements = 0;
	_num_line_elements = 0;
	_num_tri_elements = 0;
	_num_quad_elements = 0;

	model_centroid[0] = 0.0;
	model_centroid[1] = 0.0;
	model_centroid[2] = 0.0;
}

mesh_visualixer::~mesh_visualixer(){
	for (unsigned int i=0; i<_vMeshInstances.size(); i++){
		if (_vMeshInstances.at(i) == this){
			_vMeshInstances.erase(_vMeshInstances.begin() + i);
			return;
		}
	}

	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	if (_colorby != nullptr) delete[] _colorby;
	if (_color_alpha != nullptr) delete[] _color_alpha;

}

void mesh_visualixer::bind_mesh(const Mesh & mesh){

	_mesh = &mesh;

	xmax = _mesh->xmax();
	ymax = _mesh->ymax();
	zmax = _mesh->zmax();
	xmin = _mesh->xmin();
	ymin = _mesh->ymin();
	zmin = _mesh->zmin();
	//GLfloat scale = (xmax-xmin)/1.0;

	num_vertices = _mesh->nodecount();
	num_per_vertex = 7;
	num_vertex_points = 3;
	

	return;
}


void mesh_visualixer::set_test_case(){
	num_vertices = 100;
	num_per_vertex = 7;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<10; i++){
		for (unsigned int j=0; j<10; j++){
			vertices[(i*10+j)*num_per_vertex] = GLfloat(i);
			vertices[(i*10+j)*num_per_vertex + 1] = GLfloat(j);
			vertices[(i*10+j)*num_per_vertex + 2] = 10.0;
			if (j%2==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 0.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f;
			}
			else if (j%3==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 0.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f;
			}
			vertices[(i*10+j)*num_per_vertex + 6] = 1.0f;
		}
	}

	_num_point_elements = num_vertices;
	_num_line_elements = 18;
	elements = new GLuint[_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element];
	unsigned int line_element_offset = _num_point_elements*_num_per_point_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}
	// add in the line elements
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + i*_num_per_line_element] = i;
		elements[line_element_offset + i*_num_per_line_element+1] = (i+1);
	}
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + (i+9)*_num_per_line_element] = i*10;
		elements[line_element_offset + (i+9)*_num_per_line_element + 1] = (i+1)*10;
	}


	model_centroid[0] = 4.5;
	model_centroid[1] = 4.5;
	model_centroid[2] = 10.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}

void mesh_visualixer::onPrepareData(){

	MeshNode nd;
	MeshElement elem;

	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
		nd = _mesh->node(i);

		vertices[i*num_per_vertex] = nd.x();
		vertices[i*num_per_vertex + 1] = nd.y();
		vertices[i*num_per_vertex + 2] = nd.z();
		if (nd.boundary()){
			vertices[i*num_per_vertex + 3] = 1.0f;
			vertices[i*num_per_vertex + 4] = 0.0f;
			vertices[i*num_per_vertex + 5] = 0.0f;
		}
		else {
			vertices[i*num_per_vertex + 3] = 1.0f;
			vertices[i*num_per_vertex + 4] = 1.0f;
			vertices[i*num_per_vertex + 5] = 1.0f;
		}
		vertices[i*num_per_vertex + 6] = 1.0f;
	}

	// figure out how many line elements are needed
	// DYLAN_TODO: this should really be more rigorous and count for each element
	//				in the mesh
	_num_line_elements = 0;
	for (unsigned int i=0; i<_mesh->elementcount(); i++){
		elem = _mesh->element(i);
		_num_line_elements += elem.num_vertices();
	}

	_num_point_elements = num_vertices;

	elements = new GLuint[_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element];
	unsigned int line_element_offset = _num_point_elements*_num_per_point_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	
	unsigned int jp1, elements_added=0;
	for (unsigned int i=0; i<_mesh->elementcount(); i++){
		elem = _mesh->element(i);
		for (unsigned int j=0; j<elem.num_vertices(); j++){
			jp1 = (j+1)%elem.num_vertices();
			elements[line_element_offset + elements_added*_num_per_line_element] = elem.vertex_ind(j);
			elements[line_element_offset + elements_added*_num_per_line_element + 1] = elem.vertex_ind(jp1);

			elements_added++;
		}
	}


	model_centroid[0] = (xmax + xmin)/2.0;
	model_centroid[1] = (ymax + ymin)/2.0;
	model_centroid[2] = (zmax + zmin)/2.0;

}


#ifdef _TEST_
// use cmake to compile

#include "RegularMesh.hpp"

int main(int argc, char * argv[]){
	// declare vars

	// test the mesh viewer
	mesh_visualixer mymvis;
	RegularMesh mesh = RegularMesh::create_regular_grid_n(0.1, 50, 50);//, (unsigned int)30);
	mymvis.bind_mesh(mesh);
	mymvis.set_color_ramp(CRamp::DIVERGENT_9);
	mymvis.set_colorby(&mesh.x());
	//mymvis->set_test_case();
	
	mymvis.run();

	return 0;

}

#endif