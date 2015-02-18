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
	colorby = NULL;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;
	num_line_elements = 0;

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

	//delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	//if (normals != NULL) delete[] normals;
}

void mesh_visualixer::add_mesh(Static_Mesh * mesh){
	Mesh_Node nd;
	Mesh_Element elem;

	xmax = mesh->xmax();
	ymax = mesh->ymax();
	zmax = mesh->zmax();
	xmin = mesh->xmin();
	ymin = mesh->ymin();
	zmin = mesh->zmin();
	//GLfloat scale = (xmax-xmin)/1.0;

	num_vertices = mesh->nodecount();
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
		nd = mesh->node(i);

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
	}

	// figure out how many line elements are needed
	// DYLAN_TODO: this should really be more rigorous and count for each element
	//				in the mesh
	num_line_elements = 0;
	for (unsigned int i=0; i<mesh->elementcount(); i++){
		elem = mesh->element(i);
		num_line_elements += elem.num_vertices();
	}

	num_elements = num_vertices;
	num_per_element = 1;
	num_per_line_element = 2;
	elements = new GLuint[num_elements*num_per_element + num_line_elements*num_per_line_element];
	line_element_offset = num_elements*num_per_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	
	unsigned int jp1, elements_added=0;
	for (unsigned int i=0; i<mesh->elementcount(); i++){
		elem = mesh->element(i);
		for (unsigned int j=0; j<elem.num_vertices(); j++){
			jp1 = (j+1)%elem.num_vertices();
			elements[line_element_offset + elements_added*num_per_line_element] = elem.vertex_ind(j);
			elements[line_element_offset + elements_added*num_per_line_element + 1] = elem.vertex_ind(jp1);

			elements_added++;
		}
	}


	model_centroid[0] = (xmax + xmin)/2.0;
	model_centroid[1] = (ymax + ymin)/2.0;
	model_centroid[2] = (zmax + zmin)/2.0;

	return;
}


void mesh_visualixer::set_test_case(){
	num_vertices = 100;
	num_per_vertex = 6;
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
		}
	}

	num_elements = num_vertices;
	num_per_element = 1;
	num_line_elements = 18;
	num_per_line_element = 2;
	elements = new GLuint[num_elements*num_per_element + num_line_elements*num_per_line_element];
	line_element_offset = num_elements*num_per_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}
	// add in the line elements
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + i*num_per_line_element] = i;
		elements[line_element_offset + i*num_per_line_element+1] = (i+1);
	}
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + (i+9)*num_per_line_element] = i*10;
		elements[line_element_offset + (i+9)*num_per_line_element + 1] = (i+1)*10;
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

void mesh_visualixer::onRender(){
	// Create Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

	// create VBO and copy data to it
  glGenBuffers (1, &vbo);

  // visualizer specific data definitions
  // this one happens to be XYRGB
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

	// Create an element array if necessary
	if (num_elements > 0){
	    glGenBuffers(1, &ebo);
	    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (num_elements * num_per_element + num_line_elements * num_per_line_element) * sizeof(GLuint), elements, GL_STATIC_DRAW);
  }


  //cout << "total elements buffered: " << num_elements * num_per_element + num_line_elements * num_per_line_element << endl;
  // enable point size specification
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  return;
}

bool mesh_visualixer::MainLoop(){

  while(!glfwWindowShouldClose(window_ptr)){

    if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    glfwSwapBuffers(window_ptr);
    glfwPollEvents();

    // Clear the screen to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

    // Draw nodes
    glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
    // Draw lines
    glDrawElements(GL_LINES, num_line_elements*num_per_line_element , GL_UNSIGNED_INT, (void *)(line_element_offset * sizeof(GLuint)));

    //cout << "looping \r" << flush;
	}

	return 0;
}

void mesh_visualixer::onExit(){
	glDeleteProgram(shaderProgram);
  glDeleteShader(fragmentShader);
  glDeleteShader(vertexShader);

  if (num_elements > 0) glDeleteBuffers(1, &ebo);
	if (num_line_elements > 0) glDeleteBuffers(1, &lebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	glfwDestroyWindow( window_ptr );
	glfwTerminate();
	return;
}


#ifdef _TEST_
// use cmake to compile

int main(int argc, char * argv[]){
	// declare vars

	// test the mesh viewer
	mesh_visualixer * mymvis = new mesh_visualixer();
	Static_Mesh * mesh = Static_Mesh::create_regular_grid_n(0.1, 50, 50);//, (unsigned int)30);
	mymvis->add_mesh(mesh);
	mymvis->set_color_ramp(CRamp::DIVERGENT_9);
	//if (&mesh->x() == NULL) cout << "damn coloby is null" << endl;
	mymvis->set_colorby(&mesh->z());
	//mymvis->set_test_case();
	
	mymvis->run();
	delete mymvis;

	return 0;

}

#endif