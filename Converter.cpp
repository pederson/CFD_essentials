#include "Converter.hpp"

//#define _TEST_

using namespace std;

// define the functions that build a mesh from a parametric model
RegularMesh build_simple_mesh_2d(const ParametricModel2D & model,  double res, double xmin, double xmax, double ymin, double ymax, vector<double> bg_properties){

	// first create a regular grid
	RegularMesh outmesh = RegularMesh::create_regular_grid_b(res, xmin, xmax, ymin, ymax);

	// transfer background properties to the mesh
	vector<string> prop_names = model.get_phys_property_names();
	for (auto i=0; i<prop_names.size(); i++){
		outmesh.add_phys_property(prop_names.at(i), bg_properties.at(i));
	}

	vector<void *> shape_tree = model.get_object_tree();
	// perform point-in-polygon queries for each part of the parametric model
	GeometricObject2D * obj;
	for (auto i=0; i<shape_tree.size(); i++){
		obj = (GeometricObject2D *)shape_tree.at(i);
		add_shape_to_mesh(outmesh, obj, model, res);
	}

	return outmesh;
}

RegularMesh build_simple_mesh_3d(const ParametricModel3D & model, double res, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, vector<double> bg_properties){
	// first create a regular grid
	RegularMesh outmesh = RegularMesh::create_regular_grid_b(res, xmin, xmax, ymin, ymax, zmin, zmax);

	// transfer background properties to the mesh
	vector<string> prop_names = model.phys_property_names();
	for (auto i=0; i<prop_names.size(); i++){
		outmesh.add_phys_property(prop_names.at(i), bg_properties.at(i));
	}

	vector<void *> shape_tree = model.object_tree();
	// perform point-in-polygon queries for each part of the parametric model
	GeometricObject3D * obj;
	for (auto i=0; i<shape_tree.size(); i++){
		obj = (GeometricObject3D *)shape_tree.at(i);
		add_shape_to_mesh(outmesh, obj, model, res);
	}

	return outmesh;
}

void add_shape_to_mesh(RegularMesh & mesh, const GeometricObject3D * shape, const ParametricModel3D & model, double res){
	vector<string> propnames = model.phys_property_names();

	unsigned int nnodes = mesh.nodecount();
	if (shape->object_name().compare("Cylinder") == 0){

		const Cylinder * cyl = dynamic_cast<const Cylinder *>(shape);
		double rad, height;
		vertex3 center, normal;
		vector<double> shapeprops;
		const double * x, *y, *z;
		x = &mesh.x();
		y = &mesh.y();
		z = &mesh.z();

		rad = cyl->radius();
		height = cyl->height();
		center = cyl->center();
		normal = cyl->normal();
		shapeprops = cyl->phys_properties();

		double ldist, cdist;
		for (auto i=0; i<nnodes; i++){
			vertex3 pt(x[i], y[i], z[i]);
			ldist = (vertex3::cross((center-pt),normal)).norm();
			cdist = fabs(vertex3::dot((pt-center), normal));
			if (ldist <= rad && cdist <= height/2.0){
				for (unsigned int j=0; j<propnames.size(); j++){
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
				}
			}
		}

	}
	else if (shape->object_name().compare("Sphere") == 0){

		const Sphere * sph = dynamic_cast<const Sphere *>(shape);

		double rad;
		vertex3 center;
		vector<double> shapeprops;
		const double * x, *y, *z;

		x = &mesh.x();
		y = &mesh.y();
		z = &mesh.z();
		rad = sph->radius();
		center = sph->center();
		shapeprops = sph->phys_properties();


		for (auto i=0; i<nnodes; i++){
			if (vertex3::distsq(center, {x[i], y[i], z[i]}) <= rad*rad){
				for (unsigned int j=0; j<propnames.size(); j++){
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
				}
			}
		}

	}
	else if (shape->object_name().compare("Box") == 0){

		const Box * box = dynamic_cast<const Box *>(shape);

		double width, height, depth;
		vertex3 center, normal;
		vector<double> shapeprops;
		const double * x, *y, *z;

		x = &mesh.x();
		y = &mesh.y();
		z = &mesh.z();
		width = box->width();
		height = box->height();
		depth = box->depth();
		normal = box->normal();
		center = box->center();
		shapeprops = box->phys_properties();

		vertex3 sec, rel;
		for (auto i=0; i<nnodes; i++){
			vertex3 pt(x[i], y[i], z[i]);
			rel = pt-center;
			sec = rel - normal*vertex3::dot(normal, rel);
			if (sqrt(vertex3::dot(rel, normal)) <= height/2.0
				&& fabs(sec.x) <= depth/2.0 && fabs(sec.y) <= width/2.0){
				for (unsigned int j=0; j<propnames.size(); j++){
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
				}
			}
		}
		
	}
	else if (shape->object_name().compare("ParabolicDish") == 0){

		const ParabolicDish * par = dynamic_cast<const ParabolicDish *>(shape);

		double dist;
		vertex3 vertex, focus, normal, v2p;
		vector<double> shapeprops;
		const double * x, *y, *z;

		x = &mesh.x();
		y = &mesh.y();
		z = &mesh.z();
		dist = par->dist();;
		vertex = par->vertex();
		focus = par->focus();
		shapeprops = par->phys_properties();

		normal = focus-vertex;
		double d, r, fnr, v2fdist;
		v2fdist = normal.norm();
		normal.normalize();
		cout << "vertex to focus dir: " << normal.x << ", " << normal.y << ", " << normal.z << endl;
		for (auto i=0; i<nnodes; i++){
			vertex3 pt(x[i], y[i], z[i]);
			v2p = pt - vertex;
			d = vertex3::dot(v2p, normal);
			r = vertex3::cross(vertex-pt, normal).norm();
			fnr = 0.5*sqrt(d/v2fdist);
			if (r <= fnr && d <= dist && d >= 0){
				for (unsigned int j=0; j<propnames.size(); j++){
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
				}
			}
		}
		
	}
	else{
		cout << "THAT SHAPE ISNT YET IMPLEMENTED" << endl;
	}

	
	return;
}

void add_shape_to_mesh(RegularMesh & mesh, const GeometricObject2D * shape, const ParametricModel2D & model, double res){
	// convert the shape to a hull
	vector<string> propnames = model.get_phys_property_names();
	//cout << "IMM HURR" << shape.get_object_name() << endl;
	
	// do a point in polygon search for each mesh point
	unsigned int nnodes = mesh.nodecount();
	if (shape->get_object_name().compare("Gaussian2D") == 0){

		//gaussian_2d * gauss = dynamic_cast<gaussian_2d *>(shape);
		const Gaussian2D * gauss = dynamic_cast<const Gaussian2D *>(shape);
		double sx, sy, minval, amp;
		const double * x, *y;
		x = &mesh.x();
		y = &mesh.y();

		sx = gauss->sigma_x();
		sy = gauss->sigma_y();
		minval = gauss->min_val();
		amp = gauss->amplitude();
		vertex_2d cen = gauss->get_center();

		for (auto i=0; i<nnodes; i++){
			for (unsigned int j=0; j<propnames.size(); j++){
				mesh.set_phys_property(propnames.at(j), i, amp*exp(-((x[i] - cen.x)*(x[i] - cen.x)/(2*sx*sx) + (y[i] - cen.y)*(y[i] - cen.y)/(2*sy*sy))) + minval);
			}
		}

	}
	else if (shape->get_object_name().compare("Circle") == 0){

		//const Circle * circ = dynamic_cast<const Circle *>(&shape);

		const double * x, *y;
		x = &mesh.x();
		y = &mesh.y();

		map<string, double> params;
		params = shape->get_parameters();


		//rad = circ->radius();
		double rad = params.at("Radius");
		vertex_2d cen = shape->get_center();
		vector<double> shapeprops = shape->get_phys_properties();


		for (auto i=0; i<nnodes; i++){
			if ((x[i]-cen.x)*(x[i]-cen.x) + (y[i]-cen.y)*(y[i]-cen.y) <= rad*rad){
				for (unsigned int j=0; j<propnames.size(); j++){
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
				}
			}
		}

	}
	else{
		Hull shull = approximate_parametric_shape_2d(shape, res);
		MeshNode n;
		vector<double> shapeprops = shape->get_phys_properties();
		for (auto i=0; i<nnodes; i++){
			n = mesh.node(i);		// add the shape properties if it returns true
			if (shull.contains_point({n.x(), n.y()})){
				for (unsigned int j=0; j<propnames.size(); j++)
					mesh.set_phys_property(propnames.at(j), i, shapeprops.at(j));
			}
		}
	}

	
	return;
}

//Mutable_Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);

Hull approximate_parametric_shape_2d(const GeometricObject2D * shape, double res){

	vector<Point> pts_vec;
	vertex_2d cent;
	map<string, double> params;
	unsigned int numpts;

	cent = shape->get_center();
	params = shape->get_parameters();

	if (shape->get_object_name().compare("Rectangle") == 0){

			// bottom left
			pts_vec.push_back({cent.x - params.at("Width")/2.0, cent.y - params.at("Height")/2.0});
			// bottom right
			pts_vec.push_back({cent.x + params.at("Width")/2.0, cent.y - params.at("Height")/2.0});
			// top right
			pts_vec.push_back({cent.x + params.at("Width")/2.0, cent.y + params.at("Height")/2.0});
			// top left
			pts_vec.push_back({cent.x - params.at("Width")/2.0, cent.y + params.at("Height")/2.0});
	}
	else if (shape->get_object_name().compare("Circle") == 0){

			double rad = params.at("Radius");
			numpts = (unsigned int) (3.14159*2*rad/res);

			// start at right, 0 degrees, moving counter clockwise
			for (auto i=0; i<numpts; i++){
				pts_vec.push_back({cent.x + rad*cos(3.14159*2.0*(float(i)/numpts)), cent.y + rad*sin(3.14159*2.0*(float(i)/numpts))});
			}
	}
	else if(shape->get_object_name().compare("Ellipse") == 0){

			cout << "ERROR: Ellipse is not yet implemented" << endl;
			throw -1;

	}
	else if (shape->get_object_name().compare("Triangle") == 0){

			vector<vertex_2d> polyverts;
			polyverts = shape->get_vertices();

			for (auto i=0; i< polyverts.size(); i++){
				pts_vec.push_back({polyverts.at(i).x, polyverts.at(i).y});
			}
	}
	else if (shape->get_object_name().compare("Polygon") == 0){

			vector<vertex_2d> polyverts;
			polyverts = shape->get_vertices();

			for (auto i=0; i< polyverts.size(); i++){
				pts_vec.push_back({polyverts.at(i).x, polyverts.at(i).y});
			}
	}
	else{ 
		cout << "Object name not recognized" << endl;
		throw -1;
	}

	Hull approx_hull(pts_vec);
	//approx_hull->print_summary();

	return approx_hull;
}

#ifdef _TEST_

int main(int argc, char * argv[]){

	return 0;
}


#endif