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
struct vertex2{
	vertex2(){};
	vertex2(double _x, double _y){x = _x; y = _y; return;};
	double x;
	double y;

	void normalize(){
		double nval = norm();
		x = x/nval;
		y = y/nval;
	}

	double norm() const {
		return sqrt(distsq({0.0, 0.0}, *this));
	}

	vertex2 operator-(const vertex2 & v1) const{
		vertex2 vout;
		vout.x = x-v1.x;
		vout.y = y-v1.y;
		return vout;
	}

	vertex2 operator+(const vertex2 & v1) const{
		vertex2 vout;
		vout.x = x+v1.x;
		vout.y = y+v1.y;
		return vout;
	}

	vertex2 operator*(const double mult) const{
		vertex2 vout;
		vout.x = x*mult;
		vout.y = y*mult;
		return vout;
	}

	static double distsq(const vertex2 & v1, const vertex2 & v2){
		return (v1.x-v2.x)*(v1.x-v2.x) + (v1.y-v2.y)*(v1.y-v2.y);
	}

	static vertex2 orthonorm(const vertex2 & v1){
		vertex2 vout;
		vout.x = v1.y;
		vout.y = -v1.x;
		vout.normalize();
		return vout;
	}

	static double dot(const vertex2 & v1, const vertex2 & v2){
		return v1.x*v2.x + v1.y*v2.y;
	}
};

struct vertex3{
	vertex3(){};
	vertex3(double _x, double _y, double _z):x(_x), y(_y), z(_z){};
	double x;
	double y;
	double z;

	void normalize(){
		double nval = norm();
		x = x/nval;
		y = y/nval;
		z = z/nval;
	}

	double norm() const {
		return sqrt(distsq({0.0, 0.0, 0.0}, *this));
	}

	vertex3 operator-(const vertex3 & v1) const{
		vertex3 vout;
		vout.x = x-v1.x;
		vout.y = y-v1.y;
		vout.z = z-v1.z;
		return vout;
	}

	vertex3 operator+(const vertex3 & v1) const{
		vertex3 vout;
		vout.x = x+v1.x;
		vout.y = y+v1.y;
		vout.z = z+v1.z;
		return vout;
	}

	vertex3 operator*(const double mult) const{
		vertex3 vout;
		vout.x = x*mult;
		vout.y = y*mult;
		vout.z = z*mult;
		return vout;
	}

	static double distsq(const vertex3 & v1, const vertex3 & v2){
		return (v1.x-v2.x)*(v1.x-v2.x) + (v1.y-v2.y)*(v1.y-v2.y) + (v1.z-v2.z)*(v1.z-v2.z);
	}

	static vertex3 cross(const vertex3 & v1, const vertex3 & v2){
		vertex3 vout;
		vout.x = v1.y*v2.z - v1.z*v2.y;
		vout.y = v2.x*v1.z - v2.z*v1.x;
		vout.z = v1.x*v2.y - v1.y*v2.x;
		return vout;
	}

	static double dot(const vertex3 & v1, const vertex3 & v2){
		return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	}
};

//********************************** 2D GEOMETRIES ***************************
class GeometricObject2D{
public:
	GeometricObject2D();
	//~GeometricObject2D();

	// inspectors
	virtual void print_summary() const;
	std::string object_name() const {return m_object_name;};
	std::vector<double> phys_properties() const {return m_phys_properties;};
	vertex2 center() const {return m_center;};
	std::map<std::string, double> parameters() const {return m_parameters;};
	std::vector<vertex2> vertices() const {return m_vertices;};

	// mutators
	//void rotate(vertex2 point, double deg);
	virtual void translate(float delta_x, float delta_y);
	//void mirror(LineSegment blah);
	//void set_phys_property(std::string property_name, double value);

protected:

	// common data for all derived classes
	std::string m_object_name;
	vertex2 m_center;
	std::vector<double> m_phys_properties;
	//double _xmax, _xmin, _ymax, _ymin;


	// THIS DATA BELOW IS DEPRECATED AND WILL BE DELETED SOON
	// container for derived class specific data
	std::vector<std::string> m_parameter_names;
	std::map<std::string, double> m_parameters;
	// polygons only 
	std::vector<vertex2> m_vertices;
	// THIS DATA ABOVE IS DEPRECATED
};

class LineSegment : public GeometricObject2D{
public:

private:

};


class Gaussian2D : public GeometricObject2D{
public:
	Gaussian2D(double sigma_x, double sigma_y, double amplitude, double min_val, vertex2 center_);

	double sigma_x() const {return m_sigma_x;};
	double sigma_y() const {return m_sigma_y;};
	double amplitude() const {return m_amplitude;};
	double min_val() const {return m_min_val;};

	void print_summary() const;

private:
	double m_sigma_x, m_sigma_y, m_amplitude, m_min_val;


};


class Rectangle : public GeometricObject2D{
public:
	//Rectangle(double width_, double height_);
	//Rectangle(double width_, double height_, vertex2 center_);
	//Rectangle(double width_, double height_, vertex2 center_, Material mat);
	Rectangle(double width_, double height_, vertex2 center_, std::vector<double> properties);
	//~Rectangle();

	double width() const {return m_width;};
	double height() const {return m_height;};

	void print_summary() const;

private:
	double m_width, m_height;

};

class Circle : public GeometricObject2D{
public:
	//Circle(double radius_);
	//Circle(double radius_, vertex2 center_);
	//Circle(double radius_, vertex2 center_, Material mat);
	Circle(double radius_, vertex2 center_, std::vector<double> properties);
	//~Circle();

	double radius() const {return m_radius;};

	void print_summary() const;

private:
	double m_radius;

};

class Ellipse : public GeometricObject2D{
public:
	//Ellipse(double axis_major, double axis_minor, double rot_angle=0.0);
	//Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_);
	//Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, Material mat);
	Ellipse(double axis_major, double axis_minor, double rot_angle, vertex2 center_, std::vector<double> properties);
	//~Ellipse();

	double axis_major() const {return m_axis_maj;};
	double axis_minor() const {return m_axis_min;};
	double rotation_angle() const {return m_rotation_angle;};

	void print_summary() const;

private:
	double m_axis_maj, m_axis_min, m_rotation_angle;

};

class Parabola : public GeometricObject2D{
public:

	Parabola(vertex2 vertex, vertex2 focus, double dist, std::vector<double> properties);

	double dist() const {return m_dist;};
	vertex2 vertex() const {return m_vertex;};
	vertex2 focus() const {return m_focus;};

	void print_summary() const;

private:
	double m_dist;
	vertex2 m_vertex, m_focus;

};

class Triangle : public GeometricObject2D{
public:
	//Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3);
	//Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, Material mat);
	Triangle(vertex2 vert1, vertex2 vert2, vertex2 vert3, std::vector<double> properties);
	//~Triangle();

	vertex2 vert1() const {return m_v1;};
	vertex2 vert2() const {return m_v2;};
	vertex2 vert3() const {return m_v3;};


	void print_summary() const;

private:
	vertex2 m_v1, m_v2, m_v3;

};

class Polygon : public GeometricObject2D{
public:
	//Polygon(std::vector<vertex_2d> verts);
	//Polygon(std::vector<vertex_2d> verts, Material mat);
	Polygon(std::vector<vertex2> verts, std::vector<double> properties);
	//~Polygon();

	std::vector<vertex2> vertices(){return m_vertices;};

	void print_summary() const;

private:
	std::vector<vertex2> m_vertices; // ordered Clockwise

};

class RegularPolygon : public GeometricObject2D{
public:


private:
	unsigned int m_nsides;
	double m_cent_to_vert_len; // distance from center to vertex
	//alternately do center to side length

};

class ParametricModel2D : public GeometricObject2D{
public:

	ParametricModel2D();
	//~ParametricModel2D();

	void print_summary() const;

	void set_model_name(std::string mname);
	void add_physical_property(std::string property_name);
	void add_material(std::string material_name, std::vector<double> phys_props);
	void add_object(GeometricObject2D * new_object);

	void create_lattice(GeometricObject2D * new_object, vertex2 x_basis, vertex2 y_basis, unsigned int xcount, unsigned int ycount);

	//std::vector<geometric_object_2d> get_object_tree(){return ordered_object_tree;};
	std::vector<double> material(std::string material_name) const {return m_materials.at(material_name);};
	std::vector<void *> object_tree() const {return m_ordered_object_tree;};
	std::vector<std::string> phys_property_names() const {return m_phys_property_names;};

private:
	std::string m_model_name;
	//std::vector<geometric_object_2d> ordered_object_tree;
	std::vector<void *> m_ordered_object_tree;
	std::vector<std::string> m_object_tree_names;
	std::vector<std::string> m_phys_property_names;
	std::map<std::string, std::vector<double> > m_materials;

	void add_object(void * new_object, std::string object_name);

};


//********************************* 3D GEOMETRIES ****************************
class GeometricObject3D{
public:
	
	GeometricObject3D();
	//~GeometricObject3D();

	// inspectors
	virtual void print_summary() const;
	std::string object_name() const {return m_object_name;};
	std::vector<double> phys_properties() const {return m_phys_properties;};
	vertex3 center() const {return m_center;};

	// mutators
	//virtual void rotate(vertex3 point, vertex3 normal, double deg);
	//virtual void translate(float delta_x, float delta_y, float delta_z);
	//virtual void mirror(LineSegment blah);
	//void set_phys_property(std::string property_name, double value);

protected:
	// common data for all derived classes
	std::string m_object_name;
	vertex3 m_center;
	std::vector<double> m_phys_properties;
	//double m_xmax, m_xmin, m_ymax, m_ymin, m_zmax, m_zmin;
	


};

class Cylinder : public GeometricObject3D{
public:

	Cylinder(double radius, double height, vertex3 normal, vertex3 center, std::vector<double> properties);

	double radius() const {return m_radius;};
	double height() const {return m_height;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:
	double m_radius;
	double m_height;
	vertex3 m_normal;

};

class Sphere : public GeometricObject3D{
public:

	Sphere(double radius, vertex3 center, std::vector<double> properties);

	double radius() const {return m_radius;};

	void print_summary() const;

private:
	double m_radius;
};

class Ellipsoid : public GeometricObject3D{
public:

	Ellipsoid(double axis1, double axis2, double axis3, vertex3 axis1_dir, vertex3 center, std::vector<double> properties);

	double axis1() const {return m_axis1;};
	double axis2() const {return m_axis2;};
	double axis3() const {return m_axis3;};
	vertex3 axis1_dir() const {return m_axis1_dir;};

	void print_summary() const;

private:
	double m_axis1, m_axis2, m_axis3;
	vertex3 m_axis1_dir;

};

class ParabolicDish : public GeometricObject3D{
public:

	ParabolicDish(vertex3 vertex, vertex3 focus, double dist_, std::vector<double> properties);

	vertex3 vertex() const {return m_vertex;};
	vertex3 focus() const {return m_focus;};
	double dist() const {return m_dist;};

	void print_summary() const;

private:
	vertex3 m_vertex, m_focus;
	double m_dist;

};

class Box : public GeometricObject3D{
public:

	Box(double width, double height, double depth, vertex3 normal, vertex3 center, std::vector<double> properties);

	double width() const {return m_width;};
	double height() const {return m_height;};
	double depth() const {return m_depth;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:
	double m_width;
	double m_height;
	double m_depth;
	vertex3 m_normal;

};

class Prism : public GeometricObject3D{
public:

	Prism(GeometricObject2D base, double height, vertex3 normal, vertex3 center, std::vector<double> properties);

	GeometricObject2D base() const {return m_base;};
	double height() const {return m_height;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:

	GeometricObject2D m_base;
	double m_height;
	vertex3 m_normal;
};

class Cone: public GeometricObject3D{
public:

	Cone(double radius, double height, vertex3 normal, vertex3 center, std::vector<double> properties);

	double radius() const {return m_radius;};
	double height() const {return m_height;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:
	double m_height;
	double m_radius;
	vertex3 m_normal;
};

class Pyramid: public GeometricObject3D{
public:

	Pyramid(GeometricObject2D base, double height, vertex3 normal, vertex3 center, std::vector<double> properties);

	GeometricObject2D base() const {return m_base;};
	double height() const {return m_height;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:

	GeometricObject2D m_base;
	double m_height;
	vertex3 m_normal;
};

class Torus: public GeometricObject3D{
public:

	Torus(double major_radius, double minor_radius, vertex3 normal, vertex3 center, std::vector<double> properties);

	double major_radius() const {return m_major_radius;};
	double minor_radius() const {return m_minor_radius;};
	vertex3 normal() const {return m_normal;};

	void print_summary() const;

private:
	double m_major_radius;
	double m_minor_radius;
	vertex3 m_normal;
};

class TriangleMesh : public GeometricObject3D{
public:

	TriangleMesh();
	~TriangleMesh();

	void print_summary() const;
	static TriangleMesh * read_STL(std::string filename, unsigned int byte_offset=0);


	unsigned int * vertex_inds; // 3 inds per triangle
	float * vertices; //
	float * normals;
	unsigned int triangle_count, vertex_count;

private:


	//void read_STL_internal(std::string filename, unsigned int byte_offset=0);
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

	ParametricModel3D();
	//~ParametricModel3D();

	void print_summary() const;

	void set_model_name(std::string mname);
		
	void add_physical_property(std::string property_name);
	void add_material(std::string material_name, std::vector<double> phys_props);
	void add_object(GeometricObject3D * new_object);

	//void create_lattice(GeometricObject3D * new_object, vertex3 x_basis, vertex3 y_basis, vertex3 z_basis, unsigned int xcount, unsigned int ycount, unsigned int zcount);
	std::vector<double> material(std::string material_name) const {return m_materials.at(material_name);};
	std::vector<void *> object_tree() const {return m_ordered_object_tree;};
	std::vector<std::string> phys_property_names() const {return m_phys_property_names;};

	//void union();
	//void intersection();
	//void subtraction();
	// other possible 3d combinations...?

private:
	std::string m_model_name;
	std::vector<void *> m_ordered_object_tree;
	std::vector<std::string> m_object_tree_names;
	std::vector<std::string> m_phys_property_names;
	std::map<std::string, std::vector<double> > m_materials;

	void add_object(void * new_object, std::string object_name);

};

#endif
