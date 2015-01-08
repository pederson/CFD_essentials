#include "geometric_object.hpp"

#define _TEST_
using namespace std;

geometric_object_2d::geometric_object_2d(){
	object_name = "GeometricObject2D";
}

void geometric_object_2d::print_summary(){
	cout << "Shape: " << object_name << "  BASE" << endl;
}

rectangle::rectangle(float width_, float height_, vertex_2d center_, std::vector<double> properties){
	object_name = "Rectangle";
	width = width_;
	height = height_;
	center = center_;
	phys_properties = properties;
}

void rectangle::print_summary(){
	cout << "\tShape: " << object_name << " width: " << width << " height: " << height << " center: " << center.x << ", " << center.y << endl;
}

circle::circle(float radius_, vertex_2d center_, std::vector<double> properties){
	object_name = "Circle";
	radius = radius_;
	center = center_;
	phys_properties = properties;
}

void circle::print_summary(){
	cout << "\tShape: " << object_name << " radius: " << radius << " center: " << center.x << ", " << center.y << endl;
}

ellipse::ellipse(float axis_major, float axis_minor, float rot_angle, vertex_2d center_, std::vector<double> properties){
	object_name = "Ellipse";
	axis_maj = axis_major;
	axis_min = axis_minor;
	rotation_angle = rot_angle;
	center = center_;
	phys_properties = properties;
}

void ellipse::print_summary(){
	cout << "\tShape: " << object_name << " major axis: " << axis_maj << " minor axis: " << axis_min << " rotation angle: " << rotation_angle << " center: " << center.x << ", " << center.y << endl;
}

triangle::triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties){
	object_name = "Triangle";
	v1 = vert1;
	v2 = vert2;
	v3 = vert3;
	phys_properties = properties;

	center.x = (v1.x + v2.x + v3.x)/3.0;
	center.y = (v1.y + v2.y + v3.y)/3.0;
}

void triangle::print_summary(){
	cout << "\tShape: " << object_name << " vertex1: " << v1.x << ", " << v1.y << " vertex2: " << v2.x << ", " << v2.y << " vertex3: " << v3.x << ", " << v3.y << " center: " << center.x << ", " << center.y << endl;
}

polygon::polygon(std::vector<vertex_2d> verts, std::vector<double> properties){
	object_name = "Polygon";
	vertices = verts;

	center.x = 0.0;
	center.y = 0.0;
	for (auto i=0; i<vertices.size(); i++){
		center.x += vertices.at(i).x;
		center.y += vertices.at(i).y;
	}
	center.x /= vertices.size();
	center.y /= vertices.size();

	phys_properties = properties;
}

void polygon::print_summary(){
	cout << "\tShape: " << object_name << " # vertices: " << vertices.size() << endl;
}

parametric_model_2d::parametric_model_2d(){
	model_name = "DefaultModelName";
}

void parametric_model_2d::print_summary(){
	cout << "Model Name: " << model_name << endl;
	for (auto i=0; i<ordered_object_tree.size(); i++){
		ordered_object_tree.at(i).print_summary();
	}
	return;
}

void parametric_model_2d::set_model_name(std::string mname){
	model_name = mname;
	return;
}

std::vector<double> parametric_model_2d::get_material(std::string material_name){
	return materials.at(material_name);
}

void parametric_model_2d::add_physical_property(std::string property_name){
	phys_property_names.push_back(property_name);
}

void parametric_model_2d::add_material(std::string material_name, std::vector<double> phys_props){
	materials[material_name] = phys_props;
}

void parametric_model_2d::add_object(geometric_object_2d new_object){
	ordered_object_tree.push_back(new_object);
}


//************************************************************************
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

	// test 2D parametric builder
	parametric_model_2d my_param2;
	my_param2.set_model_name("Dylan Test Model");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});

	// test rectangle
	my_param2.add_object(rectangle(1.0, 2.0, vertex_2d(0.0, 0.0), my_param2.get_material("Air")));
	// test circle
	my_param2.add_object(circle(0.75, vertex_2d(0.5, 1.0), my_param2.get_material("Dielectric")));
	// test ellipse
	my_param2.add_object(ellipse(0.3, 0.2, 0.0, vertex_2d(-1.0, -1.0), my_param2.get_material("Dielectric")));
	// test triangle
	my_param2.add_object(triangle(vertex_2d(3.0, 1.0), vertex_2d(4.0, 0.0), vertex_2d(2.0, 0.0), my_param2.get_material("Air")));


	my_param2.print_summary();




	// test stl reader
	mesh_model * mytri = mesh_model::read_STL("./testfiles/brain-gear.stl");
	delete mytri;

}

#endif
