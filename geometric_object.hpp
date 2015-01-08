#ifndef GEOMETRIC_OBJECT_H
#define GEOMETRIC_OBJECT_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

#include <sys/mman.h>

enum class Material : char{
	Air, Metal, End_of_List
};

/* Make this so that I can sequentially add items on top of eachother 
*/
struct vertex_2d{
	vertex_2d(float _x, float _y){x = _x; y = _y; return;};
	float x;
	float y;
};

//********************************** 2D GEOMETRIES ***************************
class geometric_object_2d{
public:
	geometric_object_2d();
	~geometric_object_2d();

	virtual void print_summary();

	void translate(float delta_x, float delta_y);
	//void set_phys_property(std::string property_name, double value);

protected:
	std::string object_name;
	vertex_2d center;

	std::vector<double>phys_properties;

};

class rectangle : public geometric_object_2d{
public:
	rectangle(float width_, float height_);
	rectangle(float width_, float height_, vertex_2d center_);
	rectangle(float width_, float height_, vertex_2d center_, Material mat);
	rectangle(float width_, float height_, vertex_2d center_, std::vector<double> properties);
	~rectangle();

	void print_summary();

protected:
	float width, height;

};

class circle : public geometric_object_2d{
public:
	circle(float radius_);
	circle(float radius_, vertex_2d center_);
	circle(float radius_, vertex_2d center_, Material mat);
	circle(float radius_, vertex_2d center_, std::vector<double> properties);
	~circle();

	void print_summary();

protected:
	float radius;

};

class ellipse : public geometric_object_2d{
public:
	ellipse(float axis_major, float axis_minor, float rot_angle=0.0);
	ellipse(float axis_major, float axis_minor, float rot_angle, vertex_2d center_);
	ellipse(float axis_major, float axis_minor, float rot_angle, vertex_2d center_, Material mat);
	ellipse(float axis_major, float axis_minor, float rot_angle, vertex_2d center_, std::vector<double> properties);
	~ellipse();

	void print_summary();

protected:
	float axis_maj, axis_min, rotation_angle;

};

class triangle : public geometric_object_2d{
public:
	triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3);
	triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, Material mat);
	triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties);
	~triangle();

	void print_summary();

protected:
	vertex_2d v1, v2, v3;

};

class polygon : public geometric_object_2d{
public:
	polygon(std::vector<vertex_2d> verts);
	polygon(std::vector<vertex_2d> verts, Material mat);
	polygon(std::vector<vertex_2d> verts, std::vector<double> properties);
	~polygon();

	void print_summary();

protected:
	std::vector<vertex_2d> vertices; // ordered Clockwise

};

class regular_polygon : public geometric_object_2d{
public:


protected:

};

class parametric_model_2d : public geometric_object_2d{
public:

	parametric_model_2d();
	~parametric_model_2d();

	void print_summary();

	void set_model_name(std::string mname);
	void add_physical_property(std::string property_name);
	void add_object(geometric_object_2d new_object);

protected:
	std::string model_name;
	std::vector<geometric_object_2d> ordered_object_tree;
	std::vector<std::string> phys_property_names;
	std::map<std::string, std::vector<double> > materials;

};


//********************************* 3D GEOMETRIES ****************************
class geometric_object{
public:

protected:

};


class cylinder: public geometric_object{
public:

private:

};

class sphere: public geometric_object{
public:

private:
};

class box: public geometric_object{
public:

private:

};

class cone: public geometric_object{
public:

private:
};

class pyramid: public geometric_object{
public:

private:
};

class torus: public geometric_object{
public:

private:
};

class mesh_model : public geometric_object{
public:

	mesh_model();
	~mesh_model();

	static mesh_model * read_STL(char * filename, unsigned int byte_offset=0);


	unsigned int * vertex_inds; // 3 inds per triangle
	float * vertices; //
	float * normals;
	unsigned int triangle_count, vertex_count;

private:



	//float[4] model_color;

};

#pragma pack(push,1)
struct stl_tri{
	float norm_x;
	float norm_y;
	float norm_z;

	float v1_x;
	float v1_y;
	float v1_z;

	float v2_x;
	float v2_y;
	float v2_z;

	float v3_x;
	float v3_y;
	float v3_z;

	unsigned short attrib_byte_count;
};
#pragma pack(pop)

class parametric_model: public geometric_object{
public:
	//void union();
	//void intersection();
	//void subtraction();
	// other possible 3d combinations...?

protected:

private:

};

#endif
