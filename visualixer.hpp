#ifndef _VISUALIXER_H
#define _VISUALIXER_H

// this contains definitions for the openGL visualizer widget
// basically just a placeholder and reminder for now,
// but the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds


// here is the main visualizer class
class visualixer{
public:

protected:
	bool lock_rotation; // lock mouse rotations?
	float * color_ramp;

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

#endif
