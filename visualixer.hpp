#ifndef _VISUALIXER_H
#define _VISUALIXER_H

#include <iostream>
#include <vector>
//#include <thread>

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
//#include <GL/glu.h>

#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/constants.hpp>
//#include <glm/transform.hpp>

#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480
#define DEFAULT_CENTER_X 0
#define DEFAULT_CENTER_Y 0

#define VX_PI 3.14159265358979323846264338327950288

enum
{
	MOUSE_LEFT_BUTTON = 0,
	MOUSE_MIDDLE_BUTTON = 1,
	MOUSE_RIGHT_BUTTON = 2,
	MOUSE_SCROLL_UP = 3,
	MOUSE_SCROLL_DOWN = 4
};

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

	// running the visualixer
	virtual void run();


protected:
	GLFWwindow * window_ptr;
	char * window_name;
	float * color_ramp;
	float model_centroid[3]; // from [model_min, model_max]
	float screen_centroid[3]; // from [-1, 1]

	GLuint ebo, vbo, vao;
	GLuint vertexShader, fragmentShader, shaderProgram;
	GLint uniModel, uniView, uniProj;

	glm::mat4 model, view, proj;
	float rotdeg;

	// status data
	bool left_mouse_engaged; // is the left mouse button clicked?
	bool middle_mouse_engaged; // is the middle mouse button clicked?
	bool right_mouse_engaged; // is the right mouse button clicked?
	double x_upon_click, y_upon_click; // and and y positions of mouse upon click

	bool visualixer_active; // is the window currently drawn?
	bool lock_rotation; // lock mouse rotations?
	bool lock_pan; // lock mouse panning

	
	// rendering and user interaction
	//virtual void onIdle();
	virtual void onReshape(int new_width, int new_height);

	// base callbacks to interface with GLFW
	virtual void onMouseClick(int button, int action, int modifiers);
	virtual void onMouseWheel(double xoffset, double yoffset);
	virtual void onKeyboard(int key, int scancode, int action, int modifiers);
	virtual void onCursorPosition(double xpos, double ypos);

	// derived callbacks defined by me
	virtual void onMouseClickDrag(int x, int y);
	virtual void onMouseLeftDrag(int x, int y);
	virtual void onMouseRightDrag(int x, int y);
	//virtual void onMouseScrollDown(int button, int x, int y);
	//virtual void onMouseScrollUp(int button, int x, int y);
	
	virtual void onKeyDown(unsigned char key, int x, int y);
	virtual void onKeyUp(unsigned char key, int x, int y);
	//virtual void onKeyboardDown(unsigned char key, int x, int y);

	//virtual void Repaint();
	virtual void SetFullscreen(bool bFullscreen);
	virtual void Hide();
	virtual void Show();
	virtual void Close();


	// functions related to the context creation and main loop rendering
	virtual void onInit();
	virtual void onRender();
	virtual void onShaders();
	virtual void onExit();
	const GLchar * VertexShaderSource();
	const GLchar * FragmentShaderSource(); 
	bool MainLoop();

	// random shit
	virtual void run_test_triangle();

private:

	// base callbacks for GLFW functions
	static void sMouseClick(GLFWwindow * window, int button, int action, int modifiers);
	static void sMouseWheel(GLFWwindow * window, double xoffset, double yoffset);
	static void sKeyboard(GLFWwindow * window, int key, int scancode, int action, int modifiers);
	static void sCursorPosition(GLFWwindow * window, double xpos, double ypos);


	/*
	// callback functions for glut functions
	static void sClose();
	static void sReshape(int w, int h);
	static void sDisplay();
	static void sMouseClickDrag(int x, int y);
	static void sKeyUp(unsigned char key, int x, int y);
	static void sKeyDown(unsigned char key, int x, int y);
	static void sSpecialKeyUp(int key, int x, int y);
	static void sSpecialKeyDown(int key, int x, int y);
	static void sIdle();
	*/

};

//*********** here are derived classes ****************

// viewing point clouds
class cloud_visualixer : public visualixer{
public:
	cloud_visualixer();
	~cloud_visualixer();

};

// viewing a mesh
class mesh_visualixer : public visualixer{

};

// viewing a simulation
class sim_visualixer{

};

// viewing an arbitrary geometry
class geometry_visualixer{

};

/*
int init_resources(void);
void onDisplay();
void free_resources();
*/

void test_triangle(void);

#endif
