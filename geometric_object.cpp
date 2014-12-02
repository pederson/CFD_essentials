#include "geometric_object.hpp"



triangle_mesh * trangle_mesh::read_STL(char * filename, unsigned int byte_offset){
	// declare vars
	triangle_mesh * outmesh;
	int fd;
	unsigned int tricount;

	// open file and fast forward
	fd = open(filename, O_RDONLY);
	if (fd < 0){
		cout << "Error opening file in read_STL" << endl;
	throw -1;
	}
	lseek(fd, byte_offset, SEEK_SET);

	// find the triangle count
	lseek(fd, 80, SEEK_CUR); // skip the header
	read(fd, &tricount, 4); // read the triangle count


	stlmap = (char *)mmap(NULL, sizeof(stl_tri)*tricount, PROT_READ, MAP_PRIVATE, fd, 0);
	if (stlmap == MAP_FAILED){
		cout << "Failed to map stl file" << endl;
	throw -1;
	}


	// copy the triangle data into structures
	stl_tri * triangles = new stl_tri[tricount];
	memcpy(triangles, stlmap, sizeof(stl_tri)*tricount);

	// copy the structure data into the member data
	for (unsigned int i=0; i<tricount; i++){

	}

	return outmesh;
}

