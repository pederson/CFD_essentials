#include "GeometricObject.hpp"

//#define _TEST_
using namespace std;

GeometricObject2D::GeometricObject2D(){
	_object_name = "GeometricObject2D";
}

void GeometricObject2D::print_summary() const{
	//cout << "printing summary" << flush;
	cout << "\tShape: " << _object_name << "\tCenter: (" << _center.x << ", " << _center.y << ")" ;
	cout << " base print summary" << endl;
	for (auto i=0; i<parameters.size(); i++) cout << "   " << parameter_names.at(i) << ": " << parameters.at(parameter_names.at(i));
	cout << endl;
}

void GeometricObject2D::translate(float delta_x, float delta_y){
	_center.x += delta_x;
	_center.y += delta_y;
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
	cout << "\tShape: " << _object_name << "\tsigma_x: " << _sigma_x << "\tsigma_y: " << _sigma_y << "\tamplitude: " << _amplitude << "\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;

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
	cout << "\tShape: " << _object_name << "\twidth: " << _width << "\theight: " << _height << "\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;
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
	cout << "\tShape: " << _object_name << "\tradius: " << _radius << "\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;
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
	cout << "\tShape: " << _object_name << "\tmajor axis: " << _axis_maj << "\tminor axis: " << _axis_min << "\trotation angle: " << _rotation_angle << "\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;
}

Parabola::Parabola(vertex_2d vertex, vertex_2d focus, double dist, std::vector<double> properties){
	_object_name = "Parabola";
	_center	= vertex;
	phys_properties = properties;

	m_vertex = vertex;
	m_focus = focus;
	m_dist = dist;
	
}

void Parabola::print_summary() const{
	cout << "\tShape: " << _object_name << "\tvertex: (" << m_vertex.x << ", " << m_vertex.y << ")\tfocus: (" << m_focus.x << ", " << m_focus.y << ")\tdist: " << m_dist << "\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;

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
	cout << "\tShape: " << _object_name << "\tvertex1: (" << v1.x << ", " << v1.y << ")\tvertex2: (" << v2.x << ", " << v2.y << ")\tvertex3: (" << v3.x << ", " << v3.y << ")\tcenter: (" << _center.x << ", " << _center.y << ")" << endl;
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
	cout << "\tShape: " << _object_name << "\t#vertices: " << vertices.size() << endl;
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

// this leaks memory... but it's not a lot. And the current solution is an easy one that was hastily implemented
void ParametricModel2D::create_lattice(GeometricObject2D * new_object, vertex_2d basis1, vertex_2d basis2, unsigned int xcount, unsigned int ycount){
	if (new_object->get_object_name().compare("Rectangle") == 0){
		Rectangle * obj = dynamic_cast<Rectangle *> (new_object);
		Rectangle * rep;
		vertex_2d cent = obj->get_center();

		for (auto i=0; i<xcount; i++){
			for (auto j=0; j<ycount; j++){

				rep = new Rectangle(obj->width(), obj->height(), {cent.x + i*basis1.x + j*basis2.x, cent.y + i*basis1.y + j*basis2.y}, obj->get_phys_properties());
				add_object((void *) rep, "Rectangle");
			}
		}
		
	}
	else if(new_object->get_object_name().compare("Circle") == 0){
		Circle * obj = dynamic_cast<Circle *> (new_object);
		Circle * rep;
		vertex_2d cent = obj->get_center();

		for (auto i=0; i<xcount; i++){
			for (auto j=0; j<ycount; j++){

				rep = new Circle(obj->radius(), {cent.x + i*basis1.x + j*basis2.x, cent.y + i*basis1.y + j*basis2.y}, obj->get_phys_properties());
				add_object((void *) rep, "Circle");
			}
		}
		
	}
	else if(new_object->get_object_name().compare("Ellipse") == 0){
		Ellipse * obj = dynamic_cast<Ellipse *> (new_object);
		Ellipse * rep;
		vertex_2d cent = obj->get_center();

		for (auto i=0; i<xcount; i++){
			for (auto j=0; j<ycount; j++){

				rep = new Ellipse(obj->axis_major(), obj->axis_minor(), obj->rotation_angle(), {cent.x + i*basis1.x + j*basis2.x, cent.y + i*basis1.y + j*basis2.y}, obj->get_phys_properties());
				add_object((void *) rep, "Ellipse");
			}
		}
		
	}
	else if(new_object->get_object_name().compare("Triangle") == 0){
		Triangle * obj = dynamic_cast<Triangle *> (new_object);
		Triangle * rep;
		vertex_2d cent = obj->get_center();

		for (auto i=0; i<xcount; i++){
			for (auto j=0; j<ycount; j++){

				rep = new Triangle(obj->vertex1(), obj->vertex2(), obj->vertex3(), obj->get_phys_properties());
				add_object((void *) rep, "Triangle");
			}
		}
		
	}
	else if(new_object->get_object_name().compare("Polygon") == 0){
		Polygon * obj = dynamic_cast<Polygon *> (new_object);
		Polygon * rep;
		vertex_2d cent = obj->get_center();

		for (auto i=0; i<xcount; i++){
			for (auto j=0; j<ycount; j++){

				rep = new Polygon(obj->get_vertices(), obj->get_phys_properties());
				add_object((void *) rep, "Polygon");
			}
		}
		
	}
	else {
		
	}
}


//************************************************************************
GeometricObject3D::GeometricObject3D(){
	m_object_name = "GeometricObject3D";
}

void GeometricObject3D::print_summary() const{
	//cout << "printing summary" << flush;
	cout << "\tShape: " << m_object_name << "\tCenter: (" << m_center.x << ", " << m_center.y << ")" ;
	cout << " base print summary" << endl;

}

Cylinder::Cylinder(double radius, double height, vertex3 normal, vertex3 center, std::vector<double> properties){
	m_object_name = "Cylinder";
	m_radius = radius;
	m_height = height;
	m_normal = normal;
	m_center = center;
	m_phys_properties = properties;

	m_normal.normalize();
}

void Cylinder::print_summary() const{
	cout << "\tShape: " << m_object_name << "\tradius: " << m_radius << "\theight: " << m_height << "\tcenter: (" << m_center.x << ", " << m_center.y << ", " << m_center.z << ")" << endl;
}

Sphere::Sphere(double radius, vertex3 center, std::vector<double> properties){
	m_object_name = "Sphere";
	m_radius = radius;
	m_center = center;
	m_phys_properties = properties;
}

void Sphere::print_summary() const{
	cout << "\tShape: " << m_object_name << "\tradius: " << m_radius << "\tcenter: (" << m_center.x << ", " << m_center.y << ", " << m_center.z << ")" << endl;
}

Box::Box(double width, double height, double depth, vertex3 normal, vertex3 center, std::vector<double> properties){
	m_object_name = "Box";
	m_width = width;
	m_height = height;
	m_depth = depth;
	m_normal = normal;
	m_center = center;
	m_phys_properties = properties;
}

void Box::print_summary() const{
	cout << "\tShape: " << m_object_name << "\twidth: " << m_width << "\theight: " << m_height << "\tdepth: " << m_depth << "\tcenter: (" << m_center.x << ", " << m_center.y << ", " << m_center.z << ")" << endl;
}

ParabolicDish::ParabolicDish(vertex3 vertex, vertex3 focus, double dist_, std::vector<double> properties){
	m_object_name = "ParabolicDish";
	m_vertex = vertex;
	m_focus = focus;
	m_dist = dist_;
	m_phys_properties = properties;
}

void ParabolicDish::print_summary() const{
	cout << "\tShape: " << m_object_name << "\tvertex: (" << m_vertex.x << ", " << m_vertex.y << ", " << m_vertex.z << ")\tfocus " << m_focus.x << ", " << m_focus.y << ", " << m_focus.z << ")\tdist: " << m_dist << endl;
}

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

void TriangleMesh::print_summary() const{
	cout <<"\tTRIANGLE MESH" << endl;
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


ParametricModel3D::ParametricModel3D(){
	m_model_name = "DefaultModelName3D";
}

void ParametricModel3D::print_summary() const{
	cout << " " << endl;
	cout << "********* Parametric Model 3D Summary **********" << endl;
	cout << "Model Name: " << m_model_name << endl;
	for (auto i=0; i<m_ordered_object_tree.size(); i++){
		if (m_object_tree_names.at(i).compare("Cylinder") == 0){
			((Cylinder *)m_ordered_object_tree.at(i))->print_summary();
		}
		else if(m_object_tree_names.at(i).compare("Sphere") == 0){
			((Sphere *)m_ordered_object_tree.at(i))->print_summary();
		}
		else {
			((GeometricObject3D *) m_ordered_object_tree.at(i))->print_summary();
		}
	}
	cout << "*********************************************" << endl;
	cout << " " << endl;
	return;
}

void ParametricModel3D::set_model_name(std::string mname){
	m_model_name = mname;
}

void ParametricModel3D::add_physical_property(std::string property_name){
	m_phys_property_names.push_back(property_name);
}

void ParametricModel3D::add_material(std::string material_name, std::vector<double> phys_props){
	m_materials[material_name] = phys_props;
}

void ParametricModel3D::add_object(GeometricObject3D * new_object){
	if (new_object->object_name().compare("Cylinder") == 0){
		Cylinder * obj = dynamic_cast<Cylinder *> (new_object);
		add_object((void *)(obj), "Cylinder");
	}
	else if(new_object->object_name().compare("Sphere") == 0){
		Sphere * obj = dynamic_cast<Sphere *> (new_object);
		add_object((void *)obj, "Sphere");
	}
	else if(new_object->object_name().compare("Box") == 0){
		Box * obj = dynamic_cast<Box *> (new_object);
		add_object((void *)obj, "Box");
	}
	else if(new_object->object_name().compare("ParabolicDish") == 0){
		ParabolicDish * obj = dynamic_cast<ParabolicDish *> (new_object);
		add_object((void *)obj, "ParabolicDish");
	}
	else {
		add_object((void *) new_object, new_object->object_name());
	}
}

void ParametricModel3D::add_object(void * new_object, std::string object_name){
	m_object_tree_names.push_back(object_name);
	m_ordered_object_tree.push_back(new_object);
}




#ifdef _TEST_

// to compile: g++ -std=c++11 GeometricObject.cpp -o GeometricObject_test

int main(int argc, char * argv[]){
	// declare vars

	// test 2D parametric builder
	ParametricModel2D my_param2;
	my_param2.set_model_name("Dylan Test Model 2D");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});
	cout << "Finished setting materials" << endl;

	// test rectangle
	cout << "Testing rectangle..." ;
	Rectangle newrect = Rectangle(1.0, 2.0, vertex_2d(0.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newrect);
	newrect.print_summary();
	cout << "Success!" << endl;
	// test circle
	cout << "Testing circle..." ;
	Circle newcirc = Circle(0.75, vertex_2d(0.5, 1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newcirc);
	newcirc.print_summary();
	cout << "Success!" << endl;
	// test ellipse
	cout << "Testing ellipse..." ;
	Ellipse newell = Ellipse(0.3, 0.2, 0.0, vertex_2d(-1.0, -1.0), my_param2.get_material("Dielectric"));
	my_param2.add_object(&newell);
	newell.print_summary();
	cout << "Success!" << endl;
	// test triangle
	cout << "Testing triangle..." ;
	Triangle newtri = Triangle(vertex_2d(3.0, 1.0), vertex_2d(4.0, 0.0), vertex_2d(2.0, 0.0), my_param2.get_material("Air"));
	my_param2.add_object(&newtri);
	newtri.print_summary();
	cout << "Success!" << endl;
	// test gaussian 2D
	cout <<"Testing gaussian2D..." ;
	Gaussian2D newgau = Gaussian2D(0.1, 0.4, 3.5, 0.0, {0.0, 1.0});
	my_param2.add_object(&newgau);
	newgau.print_summary();
	cout << "Success!" << endl;
	// test parabola
	cout <<"Testing parabola..." ;
	Parabola newpar = Parabola({0.0, 0.2}, {0.0, 0.0}, 1.0, my_param2.get_material("Dielectric"));
	my_param2.add_object(&newpar);
	newpar.print_summary();
	cout << "Success!" << endl;
	// test polygon
	cout <<"Testing polygon..." ;
	Polygon newpol = Polygon({{3.3, 2.2}, {0.0, 2.0}, {1.0, 1.0}, {0.0, 1.0}}, my_param2.get_material("Dielectric"));
	my_param2.add_object(&newpol);
	newpol.print_summary();
	cout << "Success!" << endl;

	my_param2.print_summary();


	// test 2D parametric builder
	cout << "\n\nTesting 3D Parametric Model Builder" << endl;
	ParametricModel3D my_param3;
	my_param3.set_model_name("Dylan Test Model 3D");
	my_param3.add_physical_property("Epsilon_rel");
	my_param3.add_physical_property("Mu_rel");
	my_param3.add_material("Air", {1.0, 1.0});
	my_param3.add_material("Dielectric", {5.0, 2.0});

	// test cylinder
	cout << "Testing cylinder..." ;
	Cylinder newcyl = Cylinder(0.3, 1.0, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, my_param3.get_material("Dielectric"));
	my_param3.add_object(&newcyl);
	newcyl.print_summary();
	cout << "Success!" << endl;
	// test sphere
	cout << "Testing sphere..." ;
	Sphere newsph = Sphere(0.4, {1.0, 1.0, 1.0}, my_param3.get_material("Dielectric"));
	my_param3.add_object(&newsph);
	newsph.print_summary();
	cout << "Success!" << endl;
	// test ellipsoid
	// test parabolic dish
	// test box
	// test prism
	// test cone
	// test pyramid
	// test torus

	my_param3.print_summary();

	// test stl reader
	TriangleMesh * mytri = TriangleMesh::read_STL("./testfiles/brain-gear.stl");
	delete mytri;

}

#endif
