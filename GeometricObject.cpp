#include "GeometricObject.hpp"

//#define _TEST_
using namespace std;

geometric_object_2d::geometric_object_2d(){
	object_name = "GeometricObject2D";
}

void geometric_object_2d::print_summary() const{
	//cout << "printing summary" << flush;
	cout << "\tShape: " << object_name << "   Center: (" << center.x << ", " << center.y << ")" ;
	cout << " base print summary" << endl;
	for (auto i=0; i<parameters.size(); i++) cout << "   " << parameter_names.at(i) << ": " << parameters.at(parameter_names.at(i));
	cout << endl;
}

gaussian_2d::gaussian_2d(double sigma_x, double sigma_y, double amplitude, double min_val, vertex_2d center_){
	object_name = "Gaussian_2D";
	center = center_;
	//phys_properties = properties;

	_sigma_x = sigma_x;
	_sigma_y = sigma_y;
	_amplitude = amplitude;
	_min_val = min_val;
}

void gaussian_2d::print_summary() const{
	cout << "\tShape: " << object_name << " sigma_x: " << _sigma_x << " sigma_y: " << _sigma_y << " amplitude: " << _amplitude << " center: " << center.x << ", " << center.y << endl;

}

rectangle::rectangle(double width_, double height_, vertex_2d center_, std::vector<double> properties){
	// common parameters
	object_name = "Rectangle";
	center = center_;
	phys_properties = properties;

	// rectangle specific parameters
	parameter_names.push_back("Width");
	parameters["Width"] = width_;
	parameter_names.push_back("Height");
	parameters["Height"] = height_;
	width = width_;
	height = height_;

	
}

void rectangle::print_summary() const{
	cout << "\tShape: " << object_name << " width: " << width << " height: " << height << " center: " << center.x << ", " << center.y << endl;
}

circle::circle(double radius_, vertex_2d center_, std::vector<double> properties){
	radius = radius_;

	// common parameters
	object_name = "Circle";
	center = center_;
	phys_properties = properties;

	// circle specific parameters
	parameter_names.push_back("Radius");
	parameters["Radius"] = radius_;
}

void circle::print_summary() const{
	cout << "\tShape: " << object_name << " radius: " << radius << " center: " << center.x << ", " << center.y << endl;
}

ellipse::ellipse(double axis_major, double axis_minor, double rot_angle, vertex_2d center_, std::vector<double> properties){
	// common parameters
	object_name = "Ellipse";
	center = center_;
	phys_properties = properties;

	// ellipse specific parameters
	parameter_names.push_back("Axis_Major");
	parameters["Axis_Major"] = axis_major;
	parameter_names.push_back("Axis_Minor");
	parameters["Axis_Minor"] = axis_minor;
	parameter_names.push_back("Rotation_Angle");
	parameters["Rotation_Angle"] = rot_angle;
	axis_maj = axis_major;
	axis_min = axis_minor;
	rotation_angle = rot_angle;
}

void ellipse::print_summary() const{
	cout << "\tShape: " << object_name << " major axis: " << axis_maj << " minor axis: " << axis_min << " rotation angle: " << rotation_angle << " center: " << center.x << ", " << center.y << endl;
}

triangle::triangle(vertex_2d vert1, vertex_2d vert2, vertex_2d vert3, std::vector<double> properties){
	// common parameters
	object_name = "Triangle";
	phys_properties = properties;

	vertices.push_back(vert1);
	vertices.push_back(vert2);
	vertices.push_back(vert3);

	v1 = vert1;
	v2 = vert2;
	v3 = vert3;

	center.x = (vert1.x + vert2.x + vert3.x)/3.0;
	center.y = (vert1.y + vert2.y + vert3.y)/3.0;
}

void triangle::print_summary() const{
	cout << "\tShape: " << object_name << " vertex1: " << v1.x << ", " << v1.y << " vertex2: " << v2.x << ", " << v2.y << " vertex3: " << v3.x << ", " << v3.y << " center: " << center.x << ", " << center.y << endl;
}

polygon::polygon(std::vector<vertex_2d> verts, std::vector<double> properties){
	// common parameters
	object_name = "Polygon";

	// polygon specific parameters
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

void polygon::print_summary() const{
	cout << "\tShape: " << object_name << " # vertices: " << vertices.size() << endl;
}

parametric_model_2d::parametric_model_2d(){
	model_name = "DefaultModelName";
}

void parametric_model_2d::print_summary() const{
	cout << " " << endl;
	cout << "********* Parametric Model Summary **********" << endl;
	cout << "Model Name: " << model_name << endl;
	for (auto i=0; i<ordered_object_tree.size(); i++){
		//ordered_object_tree.at(i).print_summary();
		if (object_tree_names.at(i).compare("Rectangle") == 0){
			//geometric_object_2d * geobj = static_cast<geometric_object_2d *>(ordered_object_tree.at(i));
			//rectangle * obj = dynamic_cast<rectangle *>(geobj);
			//obj->print_summary();
			((rectangle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Circle") == 0){
			//circle * circ = (circle *) ordered_object_tree.at(i);
			//circ->print_summary();
			((circle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Ellipse") == 0){
			((ellipse *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Triangle") == 0){
			((triangle *)ordered_object_tree.at(i))->print_summary();
		}
		else if(object_tree_names.at(i).compare("Polygon") == 0){
			((polygon *) ordered_object_tree.at(i))->print_summary();
		}
		else {
			((geometric_object_2d *) ordered_object_tree.at(i))->print_summary();
		}
		//if (ordered_object_tree.at(i).get_object_name.compare("Circle")){
		//	cout << "Circle radius: " << endl;//<< ordered_object_tree.at(i).radius << endl;
		//}
	}
	cout << "*********************************************" << endl;
	cout << " " << endl;
	return;
}

void parametric_model_2d::set_model_name(std::string mname){
	model_name = mname;
	return;
}

std::vector<double> parametric_model_2d::get_material(std::string material_name) const{
	return materials.at(material_name);
}

void parametric_model_2d::add_physical_property(std::string property_name){
	phys_property_names.push_back(property_name);
}

void parametric_model_2d::add_material(std::string material_name, std::vector<double> phys_props){
	materials[material_name] = phys_props;
}

void parametric_model_2d::add_object(geometric_object_2d * new_object){
	if (new_object->get_object_name().compare("Rectangle") == 0){
		rectangle * obj = dynamic_cast<rectangle *> (new_object);
		add_object((void *)(obj), "Rectangle");
	}
	else if(new_object->get_object_name().compare("Circle") == 0){
		circle * obj = dynamic_cast<circle *> (new_object);
		add_object((void *)obj, "Circle");
	}
	else if(new_object->get_object_name().compare("Ellipse") == 0){
		ellipse * obj = dynamic_cast<ellipse *> (new_object);
		add_object((void *)obj, "Ellipse");
	}
	else if(new_object->get_object_name().compare("Triangle") == 0){
		triangle * obj = dynamic_cast<triangle *> (new_object);
		add_object((void *)obj, "Triangle");
	}
	else if(new_object->get_object_name().compare("Polygon") == 0){
		polygon * obj = dynamic_cast<polygon *> (new_object);
		add_object((void *)obj, "Polygon");
	}
	else {
		add_object((void *) new_object, new_object->get_object_name());
	}
	
}

void parametric_model_2d::add_object(void * new_object, string object_name){
	object_tree_names.push_back(object_name);
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

mesh_model * mesh_model::read_STL(string filename, unsigned int byte_offset){
	// declare vars
	mesh_model * outmesh;
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

// to compile: g++ -std=c++11 GeometricObject.cpp -o GeometricObject_test

int main(int argc, char * argv[]){
	// declare vars

	// test 2D parametric builder
	parametric_model_2d my_param2;
	my_param2.set_model_name("Dylan Test Model");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});
	cout << "Finished setting materials" << endl;

	// test rectangle
	cout << "Testing rectangle..." ;
	rectangle newrect = rectangle(1.0, 2.0, vertex_2d(0.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newrect);
	cout << "Success!" << endl;
	// test circle
	cout << "Testing circle..." ;
	circle newcirc = circle(0.75, vertex_2d(0.5, 1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newcirc);
	cout << "Success!" << endl;
	// test ellipse
	cout << "Testing ellipse..." ;
	ellipse newell = ellipse(0.3, 0.2, 0.0, vertex_2d(-1.0, -1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newell);
	cout << "Success!" << endl;
	// test triangle
	cout << "Testing triangle..." ;
	triangle newtri = triangle(vertex_2d(3.0, 1.0), vertex_2d(4.0, 0.0), vertex_2d(2.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newtri);
	cout << "Success!" << endl;

	my_param2.print_summary();




	// test stl reader
	mesh_model * mytri = mesh_model::read_STL("./testfiles/brain-gear.stl");
	delete mytri;

}

#endif
