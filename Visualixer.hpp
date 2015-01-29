#ifndef _VISUALIXER_H
#define _VISUALIXER_H

#include "ColorRamp.hpp"

#include <iostream>
#include <vector>
#include <string>

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>

#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/fast_trigonometry.hpp>
#include <glm/gtx/fast_exponential.hpp>
#include <glm/gtx/rotate_vector.hpp>
//#include <glm/transform.hpp>

#define DEFAULT_WIDTH 800
#define DEFAULT_HEIGHT 600
#define DEFAULT_CENTER_X 0
#define DEFAULT_CENTER_Y 0

#define VX_PI 3.14159265358979323846264338327950288

// this contains definitions for the openGL visualizer widget
//
// the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize parametric or mesh CAD models and move around
//  - visualize and interact with point clouds

//************** base visualixer class *********************
class visualixer{
public:

	visualixer();
	virtual ~visualixer();

	// getters and setters
	virtual std::string get_window_name() {return window_name;};
	virtual void set_window_name(std::string w_name);

	virtual void set_color_ramp(CRamp ramp_name);

	//template <class T> void set_colorby(const T * color_by);
	template <class T>
	void set_colorby(T const * color_by){
		if (colorby == NULL) colorby = new float[num_vertices];
		for (unsigned int i=0; i<num_vertices; i++) colorby[i] = float(color_by[i]);


		colorby_max = colorby[0]; colorby_min = colorby[0];
		for (auto i=1; i<num_vertices; i++){
			if (colorby[i] > colorby_max) colorby_max = colorby[i];
			if (colorby[i] < colorby_min) colorby_min = colorby[i];
		}

		return;
	}

	//virtual void set_custom_colors(rgb * colors);
	//void set_lock_rotation(bool lock_mode);

	virtual void set_test_case();

	// running the visualixer
	virtual void run();

protected:
	// model related data
	GLFWwindow * window_ptr;
	std::string window_name;
	GLfloat * vertices;
	GLuint * elements;
  	ColorRamp color_ramp;
  	float * colorby, colorby_max, colorby_min; // one per vertex
	float model_centroid[3]; // from [model_min, model_max]
	float xmin, xmax, ymin, ymax, zmin, zmax;
	unsigned int num_vertices, num_per_vertex, num_vertex_points;
	unsigned int num_elements, num_per_element; // number of triangles that may share vertices

	// shaders and buffer objects
	GLuint ebo, vbo, vao;
	GLuint vertexShader, fragmentShader, shaderProgram;
	GLint uniModel, uniView, uniProj;

	// visualizer user viewing data
	glm::mat4 model, view, proj;
	glm::vec3 eye_vec, focus_vec, up_vec;
	glm::vec3 new_eye, new_focus;
	glm::vec3 camera_side, camera_up;
	float rotdeg, zoom_level, zoom_scale, eyez_init;

	// status data
	bool left_mouse_engaged; // is the left mouse button clicked?
	bool middle_mouse_engaged; // is the middle mouse button clicked?
	bool right_mouse_engaged; // is the right mouse button clicked?
	double x_upon_click, y_upon_click; // and and y positions of mouse upon click

	bool visualixer_active; // is the window currently drawn?
	bool lock_rotation; // lock mouse rotations?
	bool lock_pan; // lock mouse panning


	// base callbacks to interface with GLFW
	virtual void onMouseClick(int button, int action, int modifiers);
	virtual void onMouseWheel(double xoffset, double yoffset);
	virtual void onKeyboard(int key, int scancode, int action, int modifiers);
	virtual void onCursorPosition(double xpos, double ypos);

	// derived callbacks defined by me
	virtual void onMouseLeftDrag(double xpos, double ypos);
	virtual void onMouseRightDrag(double xpos, double ypos);
	virtual void onMouseMiddleDrag(double xpos, double ypos);
	virtual void onKeyDown(unsigned char key, int x, int y);
	virtual void onKeyUp(unsigned char key, int x, int y);
	virtual void onReshape(int new_width, int new_height);
	virtual void SetFullscreen(bool bFullscreen);
	virtual void recalcCamera();
	virtual void cycleColorRamp();


	// functions related to the context creation and main loop rendering
	virtual const GLchar * VertexShaderSource();
	virtual const GLchar * FragmentShaderSource();
	virtual void onInit();
	virtual void onColors();
	virtual void onRender();
	virtual void onShaders();
	virtual	bool MainLoop();
	virtual void onExit();


private:

	// base callbacks for GLFW functions
	static void sMouseClick(GLFWwindow * window, int button, int action, int modifiers);
	static void sMouseWheel(GLFWwindow * window, double xoffset, double yoffset);
	static void sKeyboard(GLFWwindow * window, int key, int scancode, int action, int modifiers);
	static void sCursorPosition(GLFWwindow * window, double xpos, double ypos);

};


//************** Simuluation Visualixer *******************
class sim_visualixer{

};


#endif
