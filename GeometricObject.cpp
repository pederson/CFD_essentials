#include "GeometricObject.hpp"

//#define _TEST_
using namespace std;

GeometricObject2D::GeometricObject2D(){
	_object_name = "GeometricObject2D";
}

void GeometricObject2D::print_summary() const{
	//cout << "printing summary" << flush;
	cout << "\tShape: " << _object_name << "   Center: (" << _center.x << ", " << _center.y << ")" ;
	cout << " base print summary" << endl;
	for (auto i=0; i<parameters.size(); i++) cout << "   " << parameter_names.at(i) << ": " << parameters.at(parameter_names.at(i));
	cout << endl;
}

Gaussian2D::Gaussian2D(double sigma_x, double sigma_y, double amplitude, double min_val, vertex_2d center_){
	_object_name = "Gaussian2D";
	_center = center_;
	//phys_properties = properties;

	_sigma_x = sigma_x;
	_sigma_y = sigma_y;
	_amplitude = amplitude;
	_min_val = min_val;
}

void Gaussian2D::print_summary() const{
	cout << "\tShape: " << _object_name << " sigma_x: " << _sigma_x << " sigma_y: " << _sigma_y << " amplitude: " << _amplitude << " center: " << _center.x << ", " << _center.y << endl;

}

Rectangle::Rectangle(double width_, double height_, vertex_2d center_, std::vector<double> properties){
	// common parameters
	_object_name = "Rectangle";
	_center = center_;
	phys_properties = properties;

	// Rectangle specific parameters
	parameter_names.push_back("Width");
	parameters["Width"] = width_;
	parameter_names.push_back("Height");
	parameters["Height"] = height_;
	_width = width_;
	_height = height_;

	
}

void Rectangle::print_summary() const{
	cout << "\tShape: " << _object_name << " width: " << _width << " height: " << _height << " center: " << _center.x << ", " << _center.y << endl;
}

Circle::Circle(double radius_, vertex_2d center_, std::vector<double> properties){
	_radius = radius_;

	// common parameters
	_object_name = "Circle";
	_center = center_;
	phys_properties = properties;

	// Circle specific parameters
	parameter_names.push_back("Radius");
	parameters["Radius"] = radius_;
}

void Circle::print_summary() const{
	cout << "\tShape: " << _object_name << " radius: " << _radius << " center: " << _center.x << ", " << _center.y << endl;
}

Ellipse::Ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, std::vector<double> properties){
	// common parameters
	_object_name = "Ellipse";
	_center = center_;
	phys_properties = properties;

	// Ellipse specific parameters
	parameter_names.push_back("Axis_Major");
	parameters["Axis_Major"] = axis_major;
	parameter_names.push_back("Axis_Minor");
	parameters["Axis_Minor"] = axis_minor;
	parameter_names.push_back("Rotation_Angle");
	parameters["Rotation_Angle"] = rot_angle;
	_axis_maj = axis_major;
	_axis_min = axis_minor;
	_rotation_angle = rot_angle;
}

void Ellipse::print_summary() const{
	cout << "\tShape: " << _object_name << " major axis: " << _axis_maj << " minor axis: " << _axis_min << " rotation angle: " << _rotation_angle << " center: " << _center.x << ", " << _center.y << endl;
}

Triangle::Triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties){
	// common parameters
	_object_name = "Triangle";
	phys_properties = properties;

	vertices.push_back(vert1);
	vertices.push_back(vert2);
	vertices.push_back(vert3);

	v1 = vert1;
	v2 = vert2;
	v3 = vert3;

	_center.x = (vert1.x + vert2.x + vert3.x)/3.0;
	_center.y = (vert1.y + vert2.y + vert3.y)/3.0;
}

void Triangle::print_summary() const{
	cout << "\tShape: " << _object_name << " vertex1: " << v1.x << ", " << v1.y << " vertex2: " << v2.x << ", " << v2.y << " vertex3: " << v3.x << ", " << v3.y << " center: " << _center.x << ", " << _center.y << endl;
}

Polygon::Polygon(std::vector<vertex_2d> verts, std::vector<double> properties){
	// common parameters
	_object_name = "Polygon";

	// Polygon specific parameters
	vertices = verts;

	_center.x = 0.0;
	_center.y = 0.0;
	for (auto i=0; i<vertices.size(); i++){
		_center.x += vertices.at(i).x;
		_center.y += vertices.at(i).y;
	}
	_center.x /= vertices.size();
	_center.y /= vertices.size();

	phys_properties = properties;
}

void Polygon::print_summary() const{
	cout << "\tShape: " << _object_name << " # vertices: " << vertices.size() << endl;
}

ParametricModel2D::ParametricModel2D(){
	model_name = "DefaultModelName";
}

void ParametricModel2D::print_summary() const{
	cout << " " << endl;
	cout << "********* Parametric Model Summary **********" << endl;
	cout << "Model Name: " << model_name << endl;
	for (auto i=0; i<ordered_object_tree.size(); i++){
		//ordered_object_tree.at(i).print_summary();
		if (object_tree_names.at(i).compare("Rectangle") == 0){
			//geometric_object_2d * geobj = static_cast<geometric_object_2d *>(ordered_object_tree.at(i));
			//rectangle * obj = dynamic_cast<rectangle *>(geobj);
			//obj->print_summary();
			((Rectangle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Circle") == 0){
			//circle * circ = (circle *) ordered_object_tree.at(i);
			//circ->print_summary();
			((Circle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Ellipse") == 0){
			((Ellipse *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Triangle") == 0){
			((Triangle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Polygon") == 0){
			((Polygon *) ordered_object_tree.at(i))->print_summary();
		}
		else {
			((GeometricObject2D *) ordered_object_tree.at(i))->print_summary();
		}
		//if (ordered_object_tree.at(i).get_object_name.compare("Circle")){
		//	cout << "Circle radius: " << endl;//<< ordered_object_tree.at(i).radius << endl;
		//}
	}
	cout << "*********************************************" << endl;
	cout << " " << endl;
	return;
}

void ParametricModel2D::set_model_name(std::string mname){
	model_name = mname;
	return;
}

std::vector<double> ParametricModel2D::get_material(std::string material_name) const{
	return materials.at(material_name);
}

void ParametricModel2D::add_physical_property(std::string property_name){
	phys_property_names.push_back(property_name);
}

void ParametricModel2D::add_material(std::string material_name, std::vector<double> phys_props){
	materials[material_name] = phys_props;
}

void ParametricModel2D::add_object(GeometricObject2D * new_object){
	if (new_object->get_object_name().compare("Rectangle") == 0){
		Rectangle * obj = dynamic_cast<Rectangle *> (new_object);
		add_object((void *)(obj), "Rectangle");
	}
	else if(new_object->get_object_name().compare("Circle") == 0){
		Circle * obj = dynamic_cast<Circle *> (new_object);
		add_object((void *)obj, "Circle");
	}
	else if(new_object->get_object_name().compare("Ellipse") == 0){
		Ellipse * obj = dynamic_cast<Ellipse *> (new_object);
		add_object((void *)obj, "Ellipse");
	}
	else if(new_object->get_object_name().compare("Triangle") == 0){
		Triangle * obj = dynamic_cast<Triangle *> (new_object);
		add_object((void *)obj, "Triangle");
	}
	else if(new_object->get_object_name().compare("Polygon") == 0){
		Polygon * obj = dynamic_cast<Polygon *> (new_object);
		add_object((void *)obj, "Polygon");
	}
	else {
		add_object((void *) new_object, new_object->get_object_name());
	}
	
}

void ParametricModel2D::add_object(void * new_object, string object_name){
	object_tree_names.push_back(object_name);
	ordered_object_tree.push_back(new_object);
	
}


//************************************************************************
TriangleMesh::TriangleMesh(){
	vertices = NULL;
	normals = NULL;
	vertex_inds = NULL;
}

TriangleMesh::~TriangleMesh(){
	if(vertices != NULL) delete[] vertices;
	if(normals != NULL) delete[] normals;
	if(vertex_inds != NULL) delete[] vertex_inds;
}

TriangleMesh * TriangleMesh::read_STL(string filename, unsigned int byte_offset){
	// declare vars
	TriangleMesh * outmesh;
	int fd;
	unsigned int tricount;
	char * stlmap;

	// open file and fast forward
	fd = open(filename.c_str(), O_RDONLY);
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
	outmesh = new TriangleMesh();
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
    cout << "ruh roh! problem unmapping STL file" << endl;
    throw -1;
  }
  close(fd);

  delete[] triangles;

	return outmesh;
}



#ifdef _TEST_

// to compile: g++ -std=c++11 GeometricObject.cpp -o GeometricObject_test

int main(int argc, char * argv[]){
	// declare vars

	// test 2D parametric builder
	ParametricModel2D my_param2;
	my_param2.set_model_name("Dylan Test Model");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});
	cout << "Finished setting materials" << endl;

	// test rectangle
	cout << "Testing rectangle..." ;
	Rectangle newrect = Rectangle(1.0, 2.0, vertex_2d(0.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newrect);
	cout << "Success!" << endl;
	// test circle
	cout << "Testing circle..." ;
	Circle newcirc = Circle(0.75, vertex_2d(0.5, 1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newcirc);
	cout << "Success!" << endl;
	// test ellipse
	cout << "Testing ellipse..." ;
	Ellipse newell = Ellipse(0.3, 0.2, 0.0, vertex_2d(-1.0, -1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newell);
	cout << "Success!" << endl;
	// test triangle
	cout << "Testing triangle..." ;
	Triangle newtri = Triangle(vertex_2d(3.0, 1.0), vertex_2d(4.0, 0.0), vertex_2d(2.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newtri);
	cout << "Success!" << endl;

	my_param2.print_summary();




	// test stl reader
	TriangleMesh * mytri = TriangleMesh::read_STL("./testfiles/brain-gear.stl");
	delete mytri;

}

#endif
