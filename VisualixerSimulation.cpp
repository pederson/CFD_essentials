#include "VisualixerSimulation.hpp"

using namespace std;

vector<visualixer*> _vSimulationInstances;

simulation_visualixer::simulation_visualixer(){

	_vSimulationInstances.push_back(this);

	visualixer_active = false;
	window_name = "Simulation Visualixer";
	rotation_lock = false;
	_colorby = nullptr;
	_color_alpha = nullptr;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;

	_num_point_elements = 0;
	_num_line_elements = 0;
	_num_tri_elements = 0;
	_num_quad_elements = 0;

	_freq_Hz = 30;
	_colorby_field = "";
	_alpha_field = "";
	_cur_time_step = 0;
	_increment_val = 1;
	set_color_ramp(CRamp::DIVERGENT_9);

	model_centroid[0] = 0.0;
	model_centroid[1] = 0.0;
	model_centroid[2] = 0.0;

}

simulation_visualixer::~simulation_visualixer(){
	for (unsigned int i=0; i<_vSimulationInstances.size(); i++){
		if (_vSimulationInstances.at(i) == this){
			_vSimulationInstances.erase(_vSimulationInstances.begin() + i);
			return;
		}
	}

	//delete[] window_name;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	if (_color_alpha != nullptr) delete[] _color_alpha;
	if (_colorby != nullptr) delete[] _colorby;

}

void simulation_visualixer::bind_simulation(const SimulationData & simdata){
	_simdata = &simdata;

	_mesh = _simdata->mesh();

	xmax = _mesh->xmax();
	ymax = _mesh->ymax();
	zmax = _mesh->zmax();
	xmin = _mesh->xmin();
	ymin = _mesh->ymin();
	zmin = _mesh->zmin();

	num_vertices = _mesh->nodecount();
	num_per_vertex = 7;
	num_vertex_points = 3;
}

/*
void simulation_visualixer::set_frequency_Hz(unsigned int freq){
	_freq_Hz = freq;
}
*/

void simulation_visualixer::set_snapshot_increment(unsigned int inc){
	_increment_val = inc;
}

void simulation_visualixer::set_colorby_field(std::string fieldname){
	_colorby_field = fieldname;
}

void simulation_visualixer::set_alpha_field(std::string fieldname){
	_alpha_field = fieldname;
}

void simulation_visualixer::increment_time_step(){
	//_increment_val = 10;
	_cur_time_step = (_cur_time_step+_increment_val)%_simdata->num_time_steps();
	unsigned int ntrue = num_vertices;
	num_vertices = _mesh->nodecount();
	set_colorby(&(_simdata->get_data_at_index(_cur_time_step, _colorby_field)), false);
	if (_alpha_field.compare("") != 0) set_color_alpha(&(_simdata->get_data_at_index(_cur_time_step, _alpha_field)), false);
	num_vertices = ntrue;
	onColors();
	onAlpha();
	onRefresh();
}

void simulation_visualixer::run(){
	onInit();
	onPrepareData();
	onColors();
	onAlpha();
	onRender();
	onShaders();
	MainLoop();
	onExit();
	return;
}

void simulation_visualixer::onPrepareData(){
	
	if (_mesh->num_dims() == 1 || _mesh->num_dims() == 2){
		const Mesh * themesh = _simdata->mesh();
		MeshNode nd;
		MeshElement elem;


		vertices = new GLfloat[num_vertices*num_per_vertex];
		for (unsigned int i=0; i<num_vertices; i++){
			nd = themesh->node(i);

			vertices[i*num_per_vertex] = nd.x();
			vertices[i*num_per_vertex + 1] = nd.y();
			vertices[i*num_per_vertex + 2] = nd.z();
			if (nd.boundary()){
				vertices[i*num_per_vertex + 3] = 1.0f;
				vertices[i*num_per_vertex + 4] = 0.0f;
				vertices[i*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[i*num_per_vertex + 3] = 1.0f;
				vertices[i*num_per_vertex + 4] = 1.0f;
				vertices[i*num_per_vertex + 5] = 1.0f;
			}
			vertices[i*num_per_vertex + 6] = 1.0f;
		}

		// figure out how many line elements are needed
		// DYLAN_TODO: this should really be more rigorous and count for each element
		//				in the mesh
		_num_line_elements = 0;
		for (unsigned int i=0; i<themesh->elementcount(); i++){
			elem = themesh->element(i);
			_num_line_elements += elem.num_vertices();
		}

		_num_point_elements = num_vertices;
		elements = new GLuint[_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element];
		unsigned int line_element_offset = _num_point_elements*_num_per_point_element;
		// set the point elements
		for (unsigned int i=0; i<num_vertices; i++){
			elements[i] = i;
		}

		
		unsigned int jp1, elements_added=0;
		for (unsigned int i=0; i<themesh->elementcount(); i++){
			elem = themesh->element(i);
			for (unsigned int j=0; j<elem.num_vertices(); j++){
				jp1 = (j+1)%elem.num_vertices();
				elements[line_element_offset + elements_added*_num_per_line_element] = elem.vertex_ind(j);
				elements[line_element_offset + elements_added*_num_per_line_element + 1] = elem.vertex_ind(jp1);

				elements_added++;
			}
		}

		model_centroid[0] = (xmax + xmin)/2.0;
		model_centroid[1] = (ymax + ymin)/2.0;
		model_centroid[2] = (zmax + zmin)/2.0;

		if (_colorby_field.compare("") != 0){
			// set colorby max and min
			
			const double * dat = &(_simdata->get_data_at_index(0, _colorby_field));
			_colorby_max = dat[0]; _colorby_min = dat[0];
			for (auto t=0; t<_simdata->num_time_steps(); t++){
				dat = &(_simdata->get_data_at_index(t, _colorby_field));
				for (auto i=0; i<num_vertices; i++){
					if (dat[i] > _colorby_max) _colorby_max = dat[i];
					if (dat[i] < _colorby_min) _colorby_min = dat[i];
				}
			}
			
			set_colorby(&(_simdata->get_data_at_index(0, _colorby_field)), false);
			//cout << "set min and max colorby" << endl;
			//cout << "max: " << colorby_max << "\t min: " << colorby_min << endl;
		} 
		if (_alpha_field.compare("") != 0){
			set_color_alpha(&(_simdata->get_data_at_index(0, _alpha_field)));
		}
	}
	else if (_mesh->num_dims() == 3){
		//cout << "STARTING 3 DIMS" << endl;
		

		model_centroid[0] = (xmax + xmin)/2.0;
		model_centroid[1] = (ymax + ymin)/2.0;
		model_centroid[2] = (zmax + zmin)/2.0;

		MeshNode nd;
		MeshElement elem;

		m_xslice = model_centroid[0];
		m_yslice = model_centroid[1];
		m_zslice = model_centroid[2];
		/*
		const RegularMesh * rmesh = dynamic_cast<const RegularMesh *>(_mesh);
		m_xwidth = rmesh->res();
		m_ywidth = rmesh->res();
		m_zwidth = rmesh->res();
		*/
		// heuristic calculation
		double nx, ny, nz, rho, xvol, yvol, zvol;
		xvol = xmax-xmin;
		yvol = ymax-ymin;
		zvol = zmax-zmin;
		rho = _mesh->nodecount()/(xvol*yvol*zvol);

		m_xwidth = xvol/(rho*yvol*zvol);
		m_ywidth = yvol/(rho*xvol*zvol);
		m_zwidth = zvol/(rho*xvol*yvol);

		m_xwidth = 1.5*pow(1.0/rho, 1.0/3.0);
		m_ywidth = 1.5*m_xwidth;
		m_zwidth = 1.5*m_xwidth;
		//m_xwidth = pow(m_xwidth, 1.0/3.0);
		//m_ywidth = pow(m_ywidth, 1.0/3.0);
		//m_zwidth = pow(m_zwidth, 1.0/3.0);
		
		//cout << "CALCULATED THE WIDTHS: ";
		//cout << "xw: " << m_xwidth << endl;
		//cout << "yw: " << m_ywidth << endl;
		//cout << "zw: " << m_zwidth << endl;

		// find which points in the x, y plane
		m_subset.assign(_mesh->nodecount(), false);
		Point a1(m_xslice-m_xwidth/2, ymax);
		Point a2(m_xslice+m_xwidth/2, ymax);
		Point a3(m_xslice+m_xwidth/2, ymin);
		Point a4(m_xslice-m_xwidth/2, ymin);
		Hull hxy({a1, a2, a3, a4});
		for (auto i=0; i<_mesh->nodecount(); i++){
			nd = _mesh->node(i);
			if (hxy.contains_point(Point(nd.x(), nd.y()))) m_subset.at(i) = true;
		}
		//cout << "FINISHED XY plane" << endl;

		
		// find which points in the y, z plane
		Point b1(m_zslice-m_zwidth/2, ymax);
		Point b2(m_zslice+m_zwidth/2, ymax);
		Point b3(m_zslice+m_zwidth/2, ymin);
		Point b4(m_zslice-m_zwidth/2, ymin);
		Hull hyz({b1, b2, b3, b4});
		for (auto i=0; i<_mesh->nodecount(); i++){
			nd = _mesh->node(i);
			if (hyz.contains_point(Point(nd.z(), nd.y()))) m_subset.at(i) = true;
		}
		//cout << "FINISHED YZ plane" << endl;
		

		// find which points in the y, x plane
		Point c1(m_yslice+m_ywidth/2, xmax);
		Point c2(m_yslice-m_ywidth/2, xmax);
		Point c3(m_yslice-m_ywidth/2, xmin);
		Point c4(m_yslice+m_ywidth/2, xmin);
		Hull hyx({c1, c2, c3, c4});
		for (auto i=0; i<_mesh->nodecount(); i++){
			nd = _mesh->node(i);
			if (hyx.contains_point(Point(nd.y(), nd.x()))) m_subset.at(i) = true;
		}
		//cout << "FINISHED YX plane" << endl;
		

		// find the number of true points
		num_vertices = 0;
		for (auto i=0; i<_mesh->nodecount(); i++){
			if (m_subset.at(i)) num_vertices++;
		}
		vertices = new GLfloat[num_vertices*num_per_vertex];
		unsigned int vct=0;
		for (unsigned int i=0; i<_mesh->nodecount(); i++){
			if (!m_subset.at(i)) continue;
			nd = _mesh->node(i);

			vertices[vct*num_per_vertex] = nd.x();
			vertices[vct*num_per_vertex + 1] = nd.y();
			vertices[vct*num_per_vertex + 2] = nd.z();
			if (nd.boundary()){
				vertices[vct*num_per_vertex + 3] = 1.0f;
				vertices[vct*num_per_vertex + 4] = 0.0f;
				vertices[vct*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[vct*num_per_vertex + 3] = 1.0f;
				vertices[vct*num_per_vertex + 4] = 1.0f;
				vertices[vct*num_per_vertex + 5] = 1.0f;
			}
			vertices[vct*num_per_vertex + 6] = 1.0f;
			vct++;
		}
		//cout << "FINISHED FILLING VERTICES" << endl;

		// figure out how many line elements are needed
		// DYLAN_TODO: this should really be more rigorous and count for each element
		//				in the mesh
		_num_line_elements = 0;
		/*
		for (unsigned int i=0; i<_mesh->elementcount(); i++){
			elem = _mesh->element(i);
			_num_line_elements += elem.num_vertices();
		}
		*/

		_num_point_elements = num_vertices;

		elements = new GLuint[_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element];
		unsigned int line_element_offset = _num_point_elements*_num_per_point_element;
		// set the point elements
		for (unsigned int i=0; i<num_vertices; i++){
			elements[i] = i;
		}
		//cout << "FINISHED POINT ELEMENTS" << endl;

		/*
		unsigned int jp1, elements_added=0;
		for (unsigned int i=0; i<_mesh->elementcount(); i++){
			elem = _mesh->element(i);
			for (unsigned int j=0; j<elem.num_vertices(); j++){
				jp1 = (j+1)%elem.num_vertices();
				elements[line_element_offset + elements_added*_num_per_line_element] = elem.vertex_ind(j);
				elements[line_element_offset + elements_added*_num_per_line_element + 1] = elem.vertex_ind(jp1);

				elements_added++;
			}
		}
		*/

		//cout << "ALOMOST DONE" << endl;
		if (_colorby_field.compare("") != 0){
			// set colorby max and min
			//cout << "COLORBY FIELD: " << _colorby_field << endl;
			
			const double * dat = &(_simdata->get_data_at_index(0, _colorby_field));
			_colorby_max = dat[0]; _colorby_min = dat[0];
			for (auto t=0; t<_simdata->num_time_steps(); t++){
				dat = &(_simdata->get_data_at_index(t, _colorby_field));
				for (auto i=0; i<_mesh->nodecount(); i++){
					if (dat[i] > _colorby_max) _colorby_max = dat[i];
					if (dat[i] < _colorby_min) _colorby_min = dat[i];
				}
			}
			unsigned int ntrue = num_vertices;
			num_vertices = _mesh->nodecount();
			set_colorby(&(_simdata->get_data_at_index(0, _colorby_field)), false);
			num_vertices = ntrue;
			//cout << "set min and max colorby" << endl;
			//cout << "max: " << colorby_max << "\t min: " << colorby_min << endl;
		} 
		if (_alpha_field.compare("") != 0){
			unsigned int ntrue = num_vertices;
			num_vertices = _mesh->nodecount();
			set_color_alpha(&(_simdata->get_data_at_index(0, _alpha_field)));
			num_vertices = ntrue;
		}
		//cout << "YEP HERE" << endl;
	}

}

void simulation_visualixer::onColors(){
	if (_colorby == nullptr) return;
	if (_colorby_max - _colorby_min == 0.0) return;
	// modify the vertex array to incorporate user-defined colors


	rgb ptcolor;
	if (m_subset.size() > 0){
		unsigned int vct = 0;
		for (auto i=0; i<_mesh->nodecount(); i++){
			if (! m_subset.at(i)) continue;
			ptcolor = color_ramp.get_ramp_color(float((_colorby[i])/(_colorby_max - _colorby_min)));
			vertices[(vct+1)*num_per_vertex - 4] = ptcolor.R;
			vertices[(vct+1)*num_per_vertex - 3] = ptcolor.G;
			vertices[(vct+1)*num_per_vertex - 2] = ptcolor.B;
			vct++;
		}
	}
	else{
		for (auto i=0; i<num_vertices; i++){
			ptcolor = color_ramp.get_ramp_color(float((_colorby[i])/(_colorby_max - _colorby_min)));
			vertices[(i+1)*num_per_vertex - 4] = ptcolor.R;
			vertices[(i+1)*num_per_vertex - 3] = ptcolor.G;
			vertices[(i+1)*num_per_vertex - 2] = ptcolor.B;
		}
	}

	return;
}


void simulation_visualixer::onAlpha(){
	if (_color_alpha == nullptr) return;

	if (m_subset.size() > 0){
		unsigned int vct=0;
		for (auto i=0; i<_mesh->nodecount(); i++){
			if (! m_subset.at(i)) continue;
			vertices[(vct+1)*num_per_vertex - 1] = _color_alpha[i];
			vct++;
		}
	}
	else{
		for (auto i=0; i<num_vertices; i++){
			vertices[(i+1)*num_per_vertex - 1] = _color_alpha[i];
		}
	}

}


bool simulation_visualixer::MainLoop(){


	double lastTime = glfwGetTime(), curTime;
  while(!glfwWindowShouldClose(window_ptr)){

  	curTime = glfwGetTime();
  	if (curTime-lastTime >= 1.0/double(_freq_Hz)){ 
  		//cout << "I should be incrementing: " << curTime << "\r" << flush;
  		increment_time_step();
  		lastTime = curTime;
  	}
  	

    if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    glfwSwapBuffers(window_ptr);
    glfwPollEvents();

    // Clear the screen to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

    // Draw nodes
    glDrawElements(GL_POINTS, _num_point_elements*_num_per_point_element , GL_UNSIGNED_INT, (void *)((0) * sizeof(GLuint)));
    glDrawElements(GL_LINES, _num_line_elements*_num_per_line_element , GL_UNSIGNED_INT, (void *)((_num_point_elements*_num_per_point_element) * sizeof(GLuint)));
    glDrawElements(GL_TRIANGLES, _num_tri_elements*_num_per_tri_element , GL_UNSIGNED_INT, (void *)((_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element) * sizeof(GLuint)));
    glDrawElements(GL_QUADS, _num_quad_elements*_num_per_quad_element , GL_UNSIGNED_INT, (void *)((_num_point_elements*_num_per_point_element + _num_line_elements*_num_per_line_element + _num_tri_elements*_num_per_tri_element) * sizeof(GLuint)));
    //cout << "looping \r" << flush;
	}

	return 0;
}



