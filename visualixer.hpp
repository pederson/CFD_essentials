#ifndef _VISUALIXER_H
#define _VISUALIXER_H

#include "PointCloud.hpp"
#include "mesh_class.hpp"
#include "geometric_object.hpp"

#include <iostream>
#include <vector>
//#include <thread>

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>
//#include <GL/freeglut.h>
//#include <GL/glu.h>

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

#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480
#define DEFAULT_CENTER_X 0
#define DEFAULT_CENTER_Y 0

#define VX_PI 3.14159265358979323846264338327950288

/*
enum
{
	MOUSE_LEFT_BUTTON = 0,
	MOUSE_MIDDLE_BUTTON = 1,
	MOUSE_RIGHT_BUTTON = 2,
	MOUSE_SCROLL_UP = 3,
	MOUSE_SCROLL_DOWN = 4
};
*/

// this contains definitions for the openGL visualizer widget
//
// the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds


// here is the main visualizer class
class visualixer{
public:

	visualixer();
	virtual ~visualixer();

	// getters and setters
	virtual char * get_window_name();
	virtual void set_window_name(char * w_name);
	//float * get_color_ramp(ColorRamp vcolor);
	//void set_color_ramp(char * cramp);
	//void cycle_color_ramp();
	//void set_lock_rotation(bool lock_mode);
	virtual void set_test_case();

	// running the visualixer
	virtual void run();




protected:
	// model related data
	GLFWwindow * window_ptr;
	char * window_name;
	GLfloat * vertices;
	GLuint * elements;
	float * color_ramp;
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
	glm::vec3 eye_vec, new_eye, focus_vec, new_focus, up_vec;
	float rotdeg, zoom_level, zoom_scale;

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
	//virtual void onMouseScrollDown(int button, int x, int y);
	//virtual void onMouseScrollUp(int button, int x, int y);
	virtual void onKeyDown(unsigned char key, int x, int y);
	virtual void onKeyUp(unsigned char key, int x, int y);
	//virtual void onKeyboardDown(unsigned char key, int x, int y);
	//virtual void Repaint();
	virtual void onReshape(int new_width, int new_height);
	virtual void SetFullscreen(bool bFullscreen);


	// functions related to the context creation and main loop rendering
	virtual const GLchar * VertexShaderSource();
	virtual const GLchar * FragmentShaderSource();
	virtual void onInit();
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

//*********** here are derived classes ****************

// viewing point clouds
class cloud_visualixer : public visualixer{
public:
	cloud_visualixer();
	~cloud_visualixer();

	void add_cloud(PointCloud * cloud);
	void set_test_case();
	const GLchar * VertexShaderSource();
	bool MainLoop();

};

// viewing a mesh
class mesh_visualixer : public visualixer{
public:
	mesh_visualixer();
	~mesh_visualixer();

	void add_mesh(Mesh * mesh);
	void set_test_case();
	const GLchar * VertexShaderSource();
	void onRender();
	bool MainLoop();
	void onExit();

protected:
	GLuint * line_elements;
	GLuint lebo;
	unsigned int num_line_elements, num_per_line_element, line_element_offset;
};

// viewing a simulation
class sim_visualixer{

};

// viewing an mesh model geometry (like an STL model)
class mesh_model_visualixer : public visualixer{
public:
	mesh_model_visualixer();
	~mesh_model_visualixer();

	void add_model(mesh_model * model);
	void set_test_case();
	const GLchar * VertexShaderSource();
	void onRender();
	bool MainLoop();
	void onExit();

};

// viewing an parametric geometry (like a CAD model)
class parametric_model_visualixer : public visualixer{

	parametric_model_visualixer();
	~parametric_model_visualixer();

	//void add_model(parametric_model);
	void set_test_case();
	const GLchar * VertexShaderSource();
	void onRender();
	bool MainLoop();
	void onExit();

};

#endif
