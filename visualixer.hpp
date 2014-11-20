#ifndef _VISUALIXER_H
#define _VISUALIXER_H

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>
#include <GL/freeglut.h>

#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480
#define DEFAULT_CENTER_X 0
#define DEFAULT_CENTER_Y 0

//#define PROGRAM "glversion"

//GLuint program;
//GLint attribute_coord2d;

// this contains definitions for the openGL visualizer widget
// basically just a placeholder and reminder for now,
// but the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds


// here is the main visualizer class
class visualixer{
public:

	visualixer();
	~visualixer();

	// getters and setters
	virtual char * get_window_name();
	virtual void set_window_name(char * w_name);
	//float * get_color_ramp(ColorRamp vcolor);
	//void set_color_ramp(char * cramp);
	//void cycle_color_ramp();

	// rendering and user interaction
	//virtual void onIdle(void);
	//virtual void onRender(void);
	//virtual void onResize(int new_width, int new_height);
	virtual void onInit(void);
	virtual void onExit(void);
	//virtual void onMouseDown(int button, int x, int y);
	//virtual void onMouseUp(int button, int x, int y);
	//virtual void OnMouseMove(int x, int y);
	//virtual void onLeftMouseDrag(int x, int y);
	//virtual void onMouseWheel(int new_wheel_number, int new_direction, int x, int y);
	
	//virtual void onKeyDownl(int new_key, char cAscii);
	//virtual void onKeyUp(int new_key, char cAscii);
	//virtual void Repaint();
	virtual void SetFullscreen(bool bFullscreen);
	virtual void Hide();
	virtual void Show();
	virtual void Close();
	//virtual void onRun();
	bool MainLoop();


	void run_test_triangle();

protected:
	bool visualixer_active; // is the window currently drawn?
	bool lock_rotation; // lock mouse rotations?
	char * window_name;
	float * color_ramp;

private:
	static void sClose(void);
	static void sReshape(int w, int h);
	static void sDisplay(void);
	static void sMouse(int button, int updown, int x, int y);
	static void sMouseWheel(int wheel_number, int direction, int x, int y);
	static void sMotion(int x, int y);
	static void sKeyUp(unsigned char key, int x, int y);
	static void sKeyDown(unsigned char key, int x, int y);
	static void sIdle(void);

};

//*********** here are derived classes ****************

// viewing point clouds
class cloud_visualixer{

};

// viewing a mesh
class mesh_visualixer{

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
