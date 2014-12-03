#include "geometric_object.hpp"

#define _TEST_

using namespace std;

triangle_mesh::triangle_mesh(){
	vertices = NULL;
	normals = NULL;
	vertex_inds = NULL;
}

triangle_mesh::~triangle_mesh(){
	if(vertices != NULL) delete[] vertices;
	if(normals != NULL) delete[] normals;
	if(vertex_inds != NULL) delete[] vertex_inds;
}

triangle_mesh * triangle_mesh::read_STL(char * filename, unsigned int byte_offset){
	// declare vars
	triangle_mesh * outmesh;
	int fd;
	unsigned int tricount;
	char * stlmap;

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
	outmesh = new triangle_mesh();
	outmesh->triangle_count = tricount;
	outmesh->vertices = new float[tricount*3*3];
	outmesh->normals = new float[tricount*3];
	outmesh->vertex_inds = new unsigned int[tricount*3];
	for (unsigned int i=0; i<tricount; i++){
		//cout << "I: " << i << " \r" << flush;

		outmesh->normals[i*3] = triangles[i].norm_x;
		outmesh->normals[i*3+1] = triangles[i].norm_y;
		outmesh->normals[i*3+2] = triangles[i].norm_z;

		outmesh->vertices[i*9] = triangles[i].v1_x;
		outmesh->vertices[i*9+1] = triangles[i].v1_y;
		outmesh->vertices[i*9+2] = triangles[i].v1_z;

		outmesh->vertices[i*9+3] = triangles[i].v2_x;
		outmesh->vertices[i*9+4] = triangles[i].v2_y;
		outmesh->vertices[i*9+5] = triangles[i].v2_z;

		outmesh->vertices[i*9+6] = triangles[i].v3_x;
		outmesh->vertices[i*9+7] = triangles[i].v3_y;
		outmesh->vertices[i*9+8] = triangles[i].v3_z;

		outmesh->vertex_inds[i*3] = i*3;
		outmesh->vertex_inds[i*3+1] = i*3+1;
		outmesh->vertex_inds[i*3+2] = i*3+2;
	}

	return outmesh;
}



#ifdef _TEST_

int main(int argc, char * argv[]){
	// declare vars

	// test stl reader
	triangle_mesh * mytri = triangle_mesh::read_STL("./testfiles/brain-gear.stl");

	delete mytri;

}

#endif
