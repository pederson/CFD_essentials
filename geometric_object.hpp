#ifndef GEOMETRIC_OBJECT_H
#define GEOMETRIC_OBJECT_H

#include <fstream>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>

#include <sys/mman.h>


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

class triangle_mesh : public geometric_object{
public:

	triangle_mesh();
	~triangle_mesh();

	static triangle_mesh * read_STL(char * filename, unsigned int byte_offset=0);


	unsigned int * vertex_inds; // 3 inds per triangle 
	float * vertices; // 3*triangle_count vertices
	float * normals; 
	unsigned int triangle_count;

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

class composite_model: public geometric_object{
public:
	//void union();
	//void intersection();
	//void subtraction();
	// other possible 3d combinations...?

protected:

private:

};

#endif