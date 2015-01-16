#include "model2mesh.hpp"

//#define _TEST_

using namespace std;

// define the functions that build a mesh from a parametric model
Mesh * build_simple_mesh_2d(parametric_model_2d * model,  double res, double xmin, double xmax, double ymin, double ymax, vector<double> bg_properties){

	// first create a regular grid
	Mesh * outmesh = Mesh::create_regular_grid(res, xmin, xmax, ymin, ymax);

	// transfer property names to the mesh
	vector<string> prop_names = model->get_phys_property_names();
	for (auto i=0; i<prop_names.size(); i++){
		outmesh->add_phys_property(prop_names.at(i));
	}

	// set background props
	outmesh->set_background_properties(bg_properties);

	vector<geometric_object_2d> shape_tree = model->get_object_tree();
	// perform point-in-polygon queries for each part of the parametric model
	for (auto i=0; i<shape_tree.size(); i++){
		add_shape_to_mesh(outmesh, &shape_tree.at(i), res);
	}

	return outmesh;
}

void add_shape_to_mesh(Mesh * meshmodel, geometric_object_2d * shape, double res){
	// convert the shape to a hull
	Hull * shull = approximate_parametric_shape_2d(shape, res);
	Node * n;

	// do a point in polygon search for each mesh point
	unsigned int nnodes = meshmodel->get_num_nodes();
	for (auto i=0; i<nnodes; i++){
		n = meshmodel->get_node_ptr(meshmodel->get_node_key(i));

		// add the shape properties if it returns true
		if (shull->contains_point({n->x, n->y})){
			n->phys_properties = shape->get_phys_properties();
		}
	}

	
	return;
}

//Mesh * build_delaunay_mesh_2d(parametric_model_2d * model, double xmin, double xmax, double ymin, double ymax, double res);

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