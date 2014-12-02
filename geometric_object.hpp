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

	static triangle_mesh * read_STL(char * filename);

private:

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