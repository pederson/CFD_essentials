#include "geometric_object.hpp"

//#define _TEST_

using namespace std;

mesh_model::mesh_model(){
	vertices = NULL;
	normals = NULL;
	vertex_inds = NULL;
}

mesh_model::~mesh_model(){
	if(vertices != NULL) delete[] vertices;
	if(normals != NULL) delete[] normals;
	if(vertex_inds != NULL) delete[] vertex_inds;
}

mesh_model * mesh_model::read_STL(char * filename, unsigned int byte_offset){
	// declare vars
	mesh_model * outmesh;
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


  lseek(fd, byte_offset, SEEK_CUR); // back to the beginning
	stlmap = (char *)mmap(NULL, 84 + sizeof(stl_tri)*tricount, PROT_READ, MAP_PRIVATE, fd, 0);
	if (stlmap == MAP_FAILED){
		cout << "Failed to map stl file" << endl;
	throw -1;
	}

	// copy the triangle data into structures
	stl_tri * triangles = new stl_tri[tricount];
	memcpy(triangles, &stlmap[84], sizeof(stl_tri)*tricount);

	// copy the structure data into the member data
	outmesh = new mesh_model();
	outmesh->triangle_count = tricount;
  outmesh->vertex_count = 3*tricount;
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

  /*
  for (unsigned int i=0; i<20; i++){
    cout << " TRIANGLES PREVIEW" << endl;
    cout << "vertex: " << outmesh
  }
  */

  if (munmap(stlmap, 84 + sizeof(stl_tri)*tricount) < 0){
    cout << "ruh roh! problem unmapping LAS file" << endl;
    throw -1;
  }
  close(fd);

  delete[] triangles;

	return outmesh;
}



#ifdef _TEST_

int main(int argc, char * argv[]){
	// declare vars

	// test stl reader
	mesh_model * mytri = mesh_model::read_STL("./testfiles/brain-gear.stl");

	delete mytri;

}

#endif
