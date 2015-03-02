#ifndef _GEOMETRIC_OBJECT_H
#define _GEOMETRIC_OBJECT_H

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

/*
enum class Material : char{
	Air, Metal, End_of_List
};
*/

/* Make this so that I can sequentially add items on top of eachother 
*/
struct vertex_2d{
	vertex_2d(){};
	vertex_2d(double _x, double _y){x = _x; y = _y; return;};
	double x;
	double y;
};

//********************************** 2D GEOMETRIES ***************************
class geometric_object_2d{
public:
	geometric_object_2d();
	//~geometric_object_2d();

	virtual void print_summary() const;

	std::string get_object_name() const {return object_name;};
	std::vector<double> get_phys_properties() const {return phys_properties;};
	vertex_2d get_center() const {return center;};
	std::map<std::string, double> get_parameters() const {return parameters;};
	std::vector<vertex_2d> get_vertices() const {return vertices;};
	//void translate(float delta_x, float delta_y);
	//void set_phys_property(std::string property_name, double value);

protected:

	// common data for all derived classes
	std::string object_name;
	vertex_2d center;
	std::vector<double> phys_properties;

	// container for derived class specific data
	std::vector<std::string> parameter_names;
	std::map<std::string, double> parameters;

	// polygons only 
	std::vector<vertex_2d> vertices;
};


class gaussian_2d : public geometric_object_2d{
public:
	gaussian_2d(double sigma_x, double sigma_y, double amplitude, double min_val, vertex_2d center_);

	double sigma_x() const {return _sigma_x;};
	double sigma_y() const {return _sigma_y;};
	double amplitude() const {return _amplitude;};
	double min_val() const {return _min_val;};

	void print_summary() const;

protected:
	double _sigma_x, _sigma_y, _amplitude, _min_val;


};


class rectangle : public geometric_object_2d{
public:
	//rectangle(double width_, double height_);
	//rectangle(double width_, double height_, vertex_2d center_);
	//rectangle(double width_, double height_, vertex_2d center_, Material mat);
	rectangle(double width_, double height_, vertex_2d center_, std::vector<double> properties);
	//~rectangle();

	void print_summary() const;

protected:
	double width, height;

};

class circle : public geometric_object_2d{
public:
	//circle(double radius_);
	//circle(double radius_, vertex_2d center_);
	//circle(double radius_, vertex_2d center_, Material mat);
	circle(double radius_, vertex_2d center_, std::vector<double> properties);
	//~circle();

	void print_summary() const;

protected:
	double radius;

};

class ellipse : public geometric_object_2d{
public:
	//ellipse(double axis_major, double axis_minor, double rot_angle=0.0);
	//ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_);
	//ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, Material mat);
	ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, std::vector<double> properties);
	//~ellipse();

	void print_summary() const;

protected:
	double axis_maj, axis_min, rotation_angle;

};

class triangle : public geometric_object_2d{
public:
	//triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3);
	//triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, Material mat);
	triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties);
	//~triangle();

	void print_summary() const;

protected:
	vertex_2d v1, v2, v3;

};

class polygon : public geometric_object_2d{
public:
	//polygon(std::vector<vertex_2d> verts);
	//polygon(std::vector<vertex_2d> verts, Material mat);
	polygon(std::vector<vertex_2d> verts, std::vector<double> properties);
	//~polygon();

	void print_summary() const;

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
	//~parametric_model_2d();

	void print_summary() const;

	void set_model_name(std::string mname);
	std::vector<double> get_material(std::string material_name) const;
	
	void add_physical_property(std::string property_name);
	void add_material(std::string material_name, std::vector<double> phys_props);
	void add_object(geometric_object_2d * new_object);

	//std::vector<geometric_object_2d> get_object_tree(){return ordered_object_tree;};
	std::vector<void *> get_object_tree() const {return ordered_object_tree;};
	std::vector<std::string> get_phys_property_names() const {return phys_property_names;};

protected:
	std::string model_name;
	//std::vector<geometric_object_2d> ordered_object_tree;
	std::vector<void *> ordered_object_tree;
	std::vector<std::string> object_tree_names;
	std::vector<std::string> phys_property_names;
	std::map<std::string, std::vector<double> > materials;

	void add_object(void * new_object, std::string object_name);

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

class prism: public geometric_object{
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

	static mesh_model * read_STL(std::string filename, unsigned int byte_offset=0);


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
