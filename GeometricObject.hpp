#ifndef _GEOMETRICOBJECT_H
#define _GEOMETRICOBJECT_H

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



/* Make this so that I can sequentially add items on top of eachother 
*/
struct vertex_2d{
	vertex_2d(){};
	vertex_2d(double _x, double _y){x = _x; y = _y; return;};
	double x;
	double y;
};

//********************************** 2D GEOMETRIES ***************************
class GeometricObject2D{
public:
	GeometricObject2D();
	//~GeometricObject2D();

	// inspectors
	virtual void print_summary() const;
	std::string get_object_name() const {return _object_name;};
	std::vector<double> get_phys_properties() const {return phys_properties;};
	vertex_2d get_center() const {return _center;};
	std::map<std::string, double> get_parameters() const {return parameters;};
	std::vector<vertex_2d> get_vertices() const {return vertices;};

	// mutators
	//void rotate(vertex_2d point, double deg);
	//void translate(float delta_x, float delta_y);
	//void mirror(LineSegment blah);
	//void set_phys_property(std::string property_name, double value);

protected:

	// common data for all derived classes
	std::string _object_name;
	vertex_2d _center;
	std::vector<double> phys_properties;
	//double _xmax, _xmin, _ymax, _ymin;

	// container for derived class specific data
	std::vector<std::string> parameter_names;
	std::map<std::string, double> parameters;

	// polygons only 
	std::vector<vertex_2d> vertices;
};

class LineSegment : public GeometricObject2D{
public:

private:

};


class Gaussian2D : public GeometricObject2D{
public:
	Gaussian2D(double sigma_x, double sigma_y, double amplitude, double min_val, vertex_2d center_);

	double sigma_x() const {return _sigma_x;};
	double sigma_y() const {return _sigma_y;};
	double amplitude() const {return _amplitude;};
	double min_val() const {return _min_val;};

	void print_summary() const;

private:
	double _sigma_x, _sigma_y, _amplitude, _min_val;


};


class Rectangle : public GeometricObject2D{
public:
	//Rectangle(double width_, double height_);
	//Rectangle(double width_, double height_, vertex_2d center_);
	//Rectangle(double width_, double height_, vertex_2d center_, Material mat);
	Rectangle(double width_, double height_, vertex_2d center_, std::vector<double> properties);
	//~Rectangle();

	double width() const {return _width;};
	double height() const {return _height;};

	void print_summary() const;

private:
	double _width, _height;

};

class Circle : public GeometricObject2D{
public:
	//Circle(double radius_);
	//Circle(double radius_, vertex_2d center_);
	//Circle(double radius_, vertex_2d center_, Material mat);
	Circle(double radius_, vertex_2d center_, std::vector<double> properties);
	//~Circle();

	double radius() const {return _radius;};

	void print_summary() const;

private:
	double _radius;

};

class Ellipse : public GeometricObject2D{
public:
	//Ellipse(double axis_major, double axis_minor, double rot_angle=0.0);
	//Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_);
	//Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, Material mat);
	Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, std::vector<double> properties);
	//~Ellipse();

	double axis_major() const {return _axis_maj;};
	double axis_minor() const {return _axis_min;};
	double rotation_angle() const {return _rotation_angle;};

	void print_summary() const;

private:
	double _axis_maj, _axis_min, _rotation_angle;

};

class Triangle : public GeometricObject2D{
public:
	//Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3);
	//Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, Material mat);
	Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties);
	//~Triangle();

	void print_summary() const;

private:
	vertex_2d v1, v2, v3;

};

class Polygon : public GeometricObject2D{
public:
	//Polygon(std::vector<vertex_2d> verts);
	//Polygon(std::vector<vertex_2d> verts, Material mat);
	Polygon(std::vector<vertex_2d> verts, std::vector<double> properties);
	//~Polygon();

	void print_summary() const;

private:
	std::vector<vertex_2d> vertices; // ordered Clockwise

};

class RegularPolygon : public GeometricObject2D{
public:


private:

};

class ParametricModel2D : public GeometricObject2D{
public:

	ParametricModel2D();
	//~ParametricModel2D();

	void print_summary() const;

	void set_model_name(std::string mname);
	std::vector<double> get_material(std::string material_name) const;
	
	void add_physical_property(std::string property_name);
	void add_material(std::string material_name, std::vector<double> phys_props);
	void add_object(GeometricObject2D * new_object);

	//std::vector<geometric_object_2d> get_object_tree(){return ordered_object_tree;};
	std::vector<void *> get_object_tree() const {return ordered_object_tree;};
	std::vector<std::string> get_phys_property_names() const {return phys_property_names;};

private:
	std::string model_name;
	//std::vector<geometric_object_2d> ordered_object_tree;
	std::vector<void *> ordered_object_tree;
	std::vector<std::string> object_tree_names;
	std::vector<std::string> phys_property_names;
	std::map<std::string, std::vector<double> > materials;

	void add_object(void * new_object, std::string object_name);

};


//********************************* 3D GEOMETRIES ****************************
class GeometricObject3D{
public:

protected:

};


class Cylinder: public GeometricObject3D{
public:

private:

};

class Sphere: public GeometricObject3D{
public:

private:
};

class Box: public GeometricObject3D{
public:

private:

};

class Prism: public GeometricObject3D{
public:

private:
};

class Cone: public GeometricObject3D{
public:

private:
};

class Pyramid: public GeometricObject3D{
public:

private:
};

class Torus: public GeometricObject3D{
public:

private:
};

class TriangleMesh : public GeometricObject3D{
public:

	TriangleMesh();
	~TriangleMesh();

	static TriangleMesh * read_STL(std::string filename, unsigned int byte_offset=0);


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

class ParametricModel3D: public GeometricObject3D{
public:
	//void union();
	//void intersection();
	//void subtraction();
	// other possible 3d combinations...?

protected:

private:

};

#endif
