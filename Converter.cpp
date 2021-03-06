#include "Converter.hpp"

//#define _TEST_

using namespace std;

// define the functions that build a mesh from a parametric model
Static_Mesh * build_simple_mesh_2d(parametric_model_2d * model,  double res, double xmin, double xmax, double ymin, double ymax, vector<double> bg_properties){

	// first create a regular grid
	Static_Mesh * outmesh = Static_Mesh::create_regular_grid_b(res, xmin, xmax, ymin, ymax);

	// transfer background properties to the mesh
	vector<string> prop_names = model->get_phys_property_names();
	for (auto i=0; i<prop_names.size(); i++){
		outmesh->add_phys_property(prop_names.at(i), bg_properties.at(i));
	}

	vector<void *> shape_tree = model->get_object_tree();
	// perform point-in-polygon queries for each part of the parametric model
	for (auto i=0; i<shape_tree.size(); i++){
		add_shape_to_mesh(outmesh, (geometric_object_2d *)shape_tree.at(i), model, res);
	}

	return outmesh;
}

void add_shape_to_mesh(Static_Mesh * mesh, geometric_object_2d * shape, parametric_model_2d * model, double res){
	// convert the shape to a hull
	vector<string> propnames = model->get_phys_property_names();
	
	// do a point in polygon search for each mesh point
	unsigned int nnodes = mesh->nodecount();
	if (shape->get_object_name().compare("Gaussian_2D") == 0){

		gaussian_2d * gauss = dynamic_cast<gaussian_2d *>(shape);
		double sx, sy, minval, amp;
		const double * x, *y;
		x = &mesh->x();
		y = &mesh->y();

		sx = gauss->sigma_x();
		sy = gauss->sigma_y();
		minval = gauss->min_val();
		amp = gauss->amplitude();
		vertex_2d cen = gauss->get_center();

		for (auto i=0; i<nnodes; i++){
			for (unsigned int j=0; j<propnames.size(); j++){
				mesh->set_phys_property(propnames.at(j), i, amp*exp(-((x[i] - cen.x)*(x[i] - cen.x)/(2*sx*sx) + (y[i] - cen.y)*(y[i] - cen.y)/(2*sy*sy))) + minval);
			}
		}

	}
	else{
		Hull * shull = approximate_parametric_shape_2d(shape, res);
		Mesh_Node n;
		vector<double> shapeprops = shape->get_phys_properties();
		for (auto i=0; i<nnodes; i++){
			n = mesh->node(i);		// add the shape properties if it returns true
			if (shull->contains_point({n.x(), n.y()})){
				for (unsigned int j=0; j<propnames.size(); j++)
					mesh->set_phys_property(propnames.at(j), i, shapeprops.at(j));
			}
		}
	}

	
	return;
}

//Mutable_Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);

Hull * approximate_parametric_shape_2d(geometric_object_2d * shape, double res){

	Hull * approx_hull;
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

	approx_hull = new Hull(pts_vec);
	//approx_hull->print_summary();

	return approx_hull;
}

#ifdef _TEST_

int main(int argc, char * argv[]){

	return 0;
}


#endif