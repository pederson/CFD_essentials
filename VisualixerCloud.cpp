#include "VisualixerCloud.hpp"

using namespace std;

//#define _TEST_

vector<visualixer*> _vCloudInstances;


//*************************************** PointCloud Visualixer Class *******************************************
cloud_visualixer::cloud_visualixer(){

	_vCloudInstances.push_back(this);

	visualixer_active = false;
	window_name = "Cloud Visualixer";
	lock_rotation = false;
	colorby = NULL;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;

	model_centroid[0] = 0.0;
  model_centroid[1] = 0.0;
  model_centroid[2] = 0.0;
}

cloud_visualixer::~cloud_visualixer(){
	for (unsigned int i=0; i<_vCloudInstances.size(); i++){
		if (_vCloudInstances.at(i) == this){
			_vCloudInstances.erase(_vCloudInstances.begin() + i);
			return;
		}
	}

	//delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
}

void cloud_visualixer::set_test_case(){
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
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f;
			}
		}
	}

	num_elements = num_vertices;
	num_per_element = 1;
	elements = new GLuint[num_elements*num_per_element];
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	for (unsigned int i=0; i<num_vertices; i++){
		cout << "x: " << vertices[i*num_per_vertex] << " y: " << vertices[i*num_per_vertex +1] << " z: " << vertices[i*num_per_vertex+2] << endl;
	}

	model_centroid[0] = 4.5;
	model_centroid[1] = 4.5;
	model_centroid[2] = 0.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}

void cloud_visualixer::add_cloud(PointCloud * cloud){
	const double *x, *y, *z;
	num_vertices = cloud->size();
	num_per_vertex = 6;
	num_vertex_points = 3;
	x = cloud->x_ptr();
	y = cloud->y_ptr();
	z = cloud->z_ptr();
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
		vertices[i*num_per_vertex] = x[i];
		vertices[i*num_per_vertex + 1] = y[i];
		vertices[i*num_per_vertex + 2] = z[i];

	}

	num_elements = num_vertices;
	num_per_element = 1;
	elements = new GLuint[num_elements];
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	const rgb48 * RGB = cloud->RGB_ptr();
  if (RGB != NULL){
  	
    for (unsigned int i=0; i<num_vertices; i){
      vertices[i*num_per_vertex + 3] = RGB[i].R/65,535;
      vertices[i*num_per_vertex + 4] = RGB[i].G/65,535;
      vertices[i*num_per_vertex + 5] = RGB[i].B/65,535;
    }
  }
  else {
    //rgb ptcolor;
    for (unsigned int i=0; i<num_vertices; i++){
      //ptcolor = color_ramp.get_ramp_color((cloud->z[i]-cloud->zmin)/(cloud->zmax - cloud->zmin));
      vertices[i*num_per_vertex + 3] = 0.0;
      vertices[i*num_per_vertex + 4] = 0.0;
      vertices[i*num_per_vertex + 5] = 1.0;
    }
  }

	model_centroid[0] = (cloud->x_max() + cloud->x_min())/2.0;
	model_centroid[1] = (cloud->y_max() + cloud->y_min())/2.0;
	model_centroid[2] = (cloud->z_max() + cloud->z_min())/2.0;
	xmax = cloud->x_max();
	ymax = cloud->y_max();
	zmax = cloud->z_max();
	xmin = cloud->x_min();
	ymin = cloud->y_min();
	zmin = cloud->z_min();

	// default color by Z
	colorby = new float[num_vertices];
	for (auto i=0; i<num_vertices; i++){
		colorby[i] = z[i];
	}
	colorby_max = zmax;
	colorby_min = zmin;

	return;
}

bool cloud_visualixer::MainLoop(){

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

    // Draw points
    glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
    //cout << "looping \r" << flush;
	}

	return 0;
}


#ifdef _TEST_
// use cmake to compile

int main(int argc, char * argv[]){
	// declare vars

	// test the point cloud viewer
	cloud_visualixer * mycvis = new cloud_visualixer();
	//mycvis->set_test_case();
	PointCloud * cloud = PointCloud::read_LAS("../testfiles/ComplexSRSInfo.las");
	//PointCloud * cloud = PointCloud::read_LAS("../testfiles/xyzrgb_manuscript.las");
	//PointCloud * cloud = PointCloud::read_LAS("../testfiles/LAS12_Sample_withRGB_Quick_Terrain_Modeler.las");
	mycvis->add_cloud(cloud);
    mycvis->set_color_ramp(CRamp::DIVERGENT_1);
	mycvis->run();
	delete mycvis;

	return 0;

}

#endif