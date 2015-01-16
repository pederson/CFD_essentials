#include "./visualixer.hpp"

using namespace std;

// this contains definitions for the openGL visualizer widget
//
// the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds

vector<visualixer*> _vInstances;


//*************************************** Base Visualixer Class *******************************************
visualixer::visualixer(){

	// really should keep track of instances of visualixer
	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Visualixer";
	lock_rotation = false;
	//color_ramp = NULL;
	colorby = NULL;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;

	model_centroid[0] = 0.0;
  model_centroid[1] = 0.0;
  model_centroid[2] = 0.0;
}

visualixer::~visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}

	delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
}

void visualixer::set_window_name(char * w_name){
	window_name = w_name;
	return;
}

void visualixer::set_color_ramp(CRamp ramp_name){
  color_ramp.set_ramp(ramp_name);
}

void visualixer::set_colorby(float * color_by){
	colorby = color_by;

	colorby_max = colorby[0]; colorby_min = colorby[0];
	for (auto i=1; i<num_vertices; i++){
		if (colorby[i] > colorby_max) colorby_max = colorby[i];
		if (colorby[i] < colorby_min) colorby_min = colorby[i];
	}

	return;
}

void visualixer::run(){
	onInit();
	onColors();
	onRender();
	onShaders();
	MainLoop();
	onExit();
	return;
}

void visualixer::set_test_case(){
	num_vertices = 4;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	vertices[0] = 10.0-0.5; vertices[1] = 10.0+0.5; vertices[2] = 0.0; vertices[3] = 1.0; vertices[4] = 0.0; vertices[5] = 0.0;
	vertices[6] = 10.0+0.5; vertices[7] = 10.0+0.5; vertices[8] = 0.0; vertices[9] = 0.0; vertices[10] = 1.0; vertices[11] = 0.0;
	vertices[12] = 10.0+0.5; vertices[13] = 10.0-0.5; vertices[14] = 0.0; vertices[15] = 0.0; vertices[16] = 0.0; vertices[17] = 1.0;
	vertices[18] = 10.0-0.5; vertices[19] = 10.0-0.5; vertices[20] = 0.0; vertices[21] = 1.0; vertices[22] = 1.0; vertices[23] = 1.0;

  num_elements = 2;
  num_per_element = 3;
  elements = new GLuint[num_per_element*num_elements];
  elements[0] = 0; elements[1] = 1; elements[2] = 2;
  elements[3] = 2; elements[4] = 3; elements[5] = 0;

  model_centroid[0] = 10.0;
  model_centroid[1] = 10.0;
  model_centroid[2] = 0.0;
  xmax = 0.5;
  ymax = 0.5;
  zmax = 0.0;
  xmin = -0.5;
  ymin = -0.5;
  zmin = 0.0;

	return;
}





void visualixer::onMouseClick(int button, int action, int modifiers){
	bool state;
	if (action == GLFW_PRESS) state = true;
	else if(action == GLFW_RELEASE) state = false;
	else cout << "unknown mouse state" << endl;


	switch (button){
		case GLFW_MOUSE_BUTTON_LEFT :
			left_mouse_engaged = state;
			if (state) {
				glfwGetCursorPos(window_ptr, &x_upon_click, &y_upon_click);
				new_eye = eye_vec;
			}
			else {
				eye_vec = new_eye;
				recalcCamera();
			}

			break;
		case GLFW_MOUSE_BUTTON_MIDDLE :
			middle_mouse_engaged = state;
			//cout << "middle mouse click " << (state? "on":"off")  << endl;
			break;
		case GLFW_MOUSE_BUTTON_RIGHT :
			right_mouse_engaged = state;
			if (state) {
				glfwGetCursorPos(window_ptr, &x_upon_click, &y_upon_click);
				new_eye = eye_vec;
				new_focus = focus_vec;
			}
			else {
				eye_vec = new_eye;
				focus_vec = new_focus;
				recalcCamera();
			}
			break;
		otherwise :
			cout << "unknown mouse button" << endl;
	}
	return;
}

void visualixer::onMouseWheel(double xoffset, double yoffset){
	//cout << "MOUSE WHEEL: xoffset: " << xoffset << " yoffset: " << yoffset << endl;
	float zoom_scale_new;
	if (yoffset < 0){
		zoom_level -= 0.05;

	}
	else if (yoffset > 0){
		zoom_level += 0.05;

	}
	//cout << "zoom level: " << zoom_level << endl;
	if (zoom_level > 3.0) zoom_level = 3.0; // limit the zoom in
	if (zoom_level < 0.5) zoom_level = 0.5; // limit the zoom out
	zoom_scale_new = 1/(zoom_level);
	eye_vec = focus_vec + (eye_vec-focus_vec)*(zoom_scale_new-zoom_scale + 1.0f);
	zoom_scale = zoom_scale_new;
	view = glm::lookAt(
        eye_vec,
        focus_vec,
        up_vec
    );
}

void visualixer::onKeyboard(int key, int scancode, int action, int modifiers){
	//cout << "KEYBOARD PRESS" << endl;
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){ // resets the view
		eye_vec.x = model_centroid[0];
		eye_vec.y = model_centroid[1];
		eye_vec.z = eyez_init;
		focus_vec.x = model_centroid[0];
		focus_vec.y = model_centroid[1];
		focus_vec.z = model_centroid[2];
		up_vec.x = 0.0f;
		up_vec.y = 1.0f;
		up_vec.z = 0.0f;

		recalcCamera();

		//cout << "zmax is " << zmax << endl;
		view = glm::lookAt(
	        eye_vec,
	        focus_vec,
	        up_vec
	    	);
	}
}

void visualixer::onCursorPosition(double xpos, double ypos){
	if (left_mouse_engaged) onMouseLeftDrag(xpos, ypos);
	else if(right_mouse_engaged) onMouseRightDrag(xpos, ypos);
	else if(middle_mouse_engaged) onMouseMiddleDrag(xpos, ypos);
}

//DYLAN_TODO: do this using quaternions instead
void visualixer::onMouseLeftDrag(double xpos, double ypos){
	double new_x_pos, new_y_pos;
	int width, height;
	glm::vec3 rot_vec, in_world_ip;

	glfwGetWindowSize(window_ptr, &width, &height);
	glfwGetCursorPos(window_ptr, &new_x_pos, &new_y_pos);

	rotdeg = ((new_x_pos-x_upon_click)*(new_x_pos-x_upon_click)
					 + (new_y_pos-y_upon_click)*(new_y_pos-y_upon_click))/(width*width + height*height)*VX_PI/2;


	in_world_ip = float(new_x_pos-x_upon_click)*-camera_side + float(new_y_pos-y_upon_click)*camera_up;

	// comes from the cross product
	rot_vec = glm::cross(in_world_ip, focus_vec - eye_vec);

	new_eye = focus_vec +  glm::rotate((eye_vec-focus_vec), rotdeg, rot_vec);
	view = glm::lookAt(new_eye, focus_vec, up_vec);
}

void visualixer::onMouseRightDrag(double xpos, double ypos){
	double new_x_pos, new_y_pos;
	int width, height;
	glm::vec3 in_plane, in_world_ip;
	glm::vec4 world_in_plane;
	GLfloat pan_scale, in_plane_norm;

	glfwGetWindowSize(window_ptr, &width, &height);
	glfwGetCursorPos(window_ptr, &new_x_pos, &new_y_pos);


	pan_scale = ((new_x_pos-x_upon_click)*(new_x_pos-x_upon_click) + (new_y_pos-y_upon_click)*(new_y_pos-y_upon_click))/(width*width + height*height)*((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin));
	//cout << "pan scale: " << pan_scale << endl;

	in_world_ip = float(new_x_pos-x_upon_click)*-camera_side + float(new_y_pos-y_upon_click)*camera_up;
	in_world_ip = glm::normalize(in_world_ip);

	//cout << "about to translate " << endl;
	new_eye = eye_vec + pan_scale*in_world_ip;
	new_focus = focus_vec + pan_scale*in_world_ip;
	view = glm::lookAt(new_eye, new_focus, up_vec);
	//cout << "finished translating" << endl;
}

void visualixer::onMouseMiddleDrag(double xpos, double ypos){
	cout << "dragging middle mouse button" << endl;
}

void visualixer::onKeyDown(unsigned char key, int x, int y){
	cout << "Pressed Key with unsigned char : " << key << endl;
	return;
}

void visualixer::onKeyUp(unsigned char key, int x, int y){
	cout << "Released Key with unsigned char: " << key << endl;
	return;
}

void visualixer::onReshape(int new_width, int new_height){
	cout << "reshaping" << endl;
}

void visualixer::SetFullscreen(bool bFullscreen){

}

void visualixer::recalcCamera(){
	camera_side = glm::normalize(glm::cross(focus_vec - eye_vec, up_vec - eye_vec + focus_vec));
	camera_up = glm::normalize(glm::cross(camera_side, focus_vec - eye_vec));

	/*
	cout << "focus_vec: " << focus_vec.x << ", " << focus_vec.y << ", " << focus_vec.z << endl;
	cout << "eye_vec: " << eye_vec.x << ", " << eye_vec.y << ", " << eye_vec.z << endl;
	cout << "up_vec: " << up_vec.x << ", " << up_vec.y << ", " << up_vec.z << endl;
	cout << "camera_side: " << camera_side.x << ", " << camera_side.y << ", " << camera_side.z << endl;
	cout << "camera_up: " << camera_up.x << ", " << camera_up.y << ", " << camera_up.z << endl;
	*/
}





const GLchar * visualixer::VertexShaderSource(){
	// vertex shader sources
	const GLchar* vertexSource =
	    "#version 140\n"
	    "in vec3 position;"
	    "in vec3 color;"
	    "out vec3 Color;"
	    "uniform mat4 model;"
	    "uniform mat4 view;"
    	"uniform mat4 proj;"
	    "void main() {"
	    "   Color = color;"
      "   gl_PointSize = 2.0;"
	    "   gl_Position = proj*view*model*vec4(position, 1.0);"
	    "}";
	return vertexSource;

}

const GLchar * visualixer::FragmentShaderSource(){
  // fragment shader source
	const GLchar* fragmentSource =
    "#version 140\n"
    "in vec3 Color;"
    "out vec4 outColor;"
    "void main() {"
    "   outColor = vec4(Color, 1.0);"
    "}";
  return fragmentSource;
}

void visualixer::onInit(){

	if (!glfwInit()){
		cout << "glfw Init failure " << endl;
		throw -1;
	}
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1); // This must be compatible with the installed version of OpenGL
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); THIS CAUSES A FAILURE
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);


	window_ptr = glfwCreateWindow(DEFAULT_WIDTH, DEFAULT_HEIGHT, window_name, NULL, NULL); // Windowed
	//GLFWwindow* window = glfwCreateWindow(400, 300, "OpenGL", glfwGetPrimaryMonitor(), NULL); // Fullscreen
	if ( !window_ptr ) {
		cout << "failed to create window" << endl;
        glfwTerminate();
        exit ( EXIT_FAILURE );
  }
	glfwMakeContextCurrent(window_ptr);

	glewExperimental = GL_TRUE;
	GLenum glew_status = glewInit();
    if ( GLEW_OK != glew_status ) {
        exit ( EXIT_FAILURE );
  }

  // set callbacks
  glfwSetMouseButtonCallback(window_ptr, sMouseClick);
  glfwSetScrollCallback(window_ptr, sMouseWheel);
  glfwSetKeyCallback(window_ptr, sKeyboard);
  glfwSetCursorPosCallback(window_ptr, sCursorPosition);
}

void visualixer::onColors(){
	if (colorby == NULL) return;
	if (colorby_max - colorby_min == 0.0) return;
	// modify the vertex array to incorporate user-defined colors
	rgb ptcolor;

	for (unsigned int i=0; i<num_vertices; i++){
		ptcolor = color_ramp.get_ramp_color((colorby[i]-colorby_min)/(colorby_max - colorby_min));
		vertices[i*num_per_vertex + 3] = ptcolor.R;
		vertices[i*num_per_vertex + 4] = ptcolor.G;
		vertices[i*num_per_vertex + 5] = ptcolor.B;
	}

	return;
}

void visualixer::onRender(){
	// Create Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

	// create VBO and copy data to it
  glGenBuffers (1, &vbo);

  // visualizer specific data definitions
  // this one happens to be XYRGB
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

	// Create an element array if necessary
	if (num_elements > 0){
	    glGenBuffers(1, &ebo);
	    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	    glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_elements * num_per_element * sizeof(GLuint), elements, GL_STATIC_DRAW);
  }

  // enable point size specification
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  return;
}

void visualixer::onShaders(){

  // create and compile the vertex shader
  vertexShader = glCreateShader(GL_VERTEX_SHADER);
  const GLchar * vssource = this->VertexShaderSource();
  glShaderSource(vertexShader, 1, &(vssource), NULL);
  glCompileShader(vertexShader);

  GLint status;
  GLint logLength;
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE){
		cout << "vertex shader failed to compile" << endl;
		glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH , &logLength);
		if (logLength > 1) {
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(vertexShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}

  // create and compile the fragment shader
  fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  const GLchar * fssource = this->FragmentShaderSource();
  glShaderSource(fragmentShader, 1, &(fssource), NULL);
  glCompileShader(fragmentShader);

  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE){
		cout << "fragment shader failed to compile" << endl;
		glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH , &logLength);
		if (logLength > 1) {
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(fragmentShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}


  // link the vertex and fragment shader into a shader program
  shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glBindFragDataLocation(shaderProgram, 0, "outColor");
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  // Specify the layout of the vertex data
  //cout << "num_vertices: " << num_vertices << " num_vertex_points: " << num_vertex_points << " num_per_vertex: " << num_per_vertex << endl;
  GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
  glEnableVertexAttribArray(posAttrib);
  glVertexAttribPointer(posAttrib, num_vertex_points, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat), 0);

  GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
  glEnableVertexAttribArray(colAttrib);
  glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat), (void*)(num_vertex_points * sizeof(GLfloat)));

  rotdeg = 0;
	model = glm::rotate(model, rotdeg, glm::vec3(0.0f, 0.0f, 1.0f)); // angle in radians to suppress some output
	uniModel = glGetUniformLocation(shaderProgram, "model");
	glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

	// Set up view matrix
	zoom_level = 1.0f;
	zoom_scale = 1.0f;
	up_vec = glm::vec3(0.0f, 1.0f, 0.0f);
	// calculate the eye z position so that it can view the whole scene
	if (xmax - xmin > ymax - ymin) eyez_init = (xmax-xmin) + zmax;
	else eyez_init = (ymax - ymin) + zmax;
	//cout << "centroid is: " << model_centroid[0] << ", " << model_centroid[1] << ", " << model_centroid[2] << endl;
	focus_vec = glm::vec3(model_centroid[0], model_centroid[1], model_centroid[2]);
	eye_vec = glm::vec3(model_centroid[0], model_centroid[1], eyez_init);
  view = glm::lookAt(
    eye_vec, // camera position
    focus_vec, // the position to be looking at
    up_vec  // the up vector
  );
  recalcCamera();

  uniView = glGetUniformLocation(shaderProgram, "view");
  glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

  // set up projection matrix
  proj = glm::perspective(0.785f, float(DEFAULT_WIDTH)/float(DEFAULT_HEIGHT), 0.05f, 100000.0f);
  uniProj = glGetUniformLocation(shaderProgram, "proj");
  glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

}


bool visualixer::MainLoop(){

  while(!glfwWindowShouldClose(window_ptr)){

    if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    glfwSwapBuffers(window_ptr);
    glfwPollEvents();

    // Clear the screen to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

    // Draw a triangle from the 3 vertices
    glDrawElements(GL_TRIANGLES, num_per_element*num_elements, GL_UNSIGNED_INT, 0);
    //cout << "looping \r" << flush;
	}

	return 0;
}

void visualixer::onExit(){
	glDeleteProgram(shaderProgram);
  glDeleteShader(fragmentShader);
  glDeleteShader(vertexShader);

  if (num_elements > 0) glDeleteBuffers(1, &ebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	glfwDestroyWindow( window_ptr );
	glfwTerminate();
	return;
}







void visualixer::sMouseClick(GLFWwindow * window, int button, int action, int modifiers){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == window){
			_vInstances.at(i)->onMouseClick(button, action, modifiers);
			return;
		}
	}
	return;
}


void visualixer::sMouseWheel(GLFWwindow * window, double xoffset, double yoffset){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == window){
			_vInstances.at(i)->onMouseWheel(xoffset, yoffset);
			return;
		}
	}
	return;
}

void visualixer::sKeyboard(GLFWwindow * window, int key, int scancode, int action, int modifiers){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == window){
			_vInstances.at(i)->onKeyboard(key, scancode, action, modifiers);
			return;
		}
	}
	return;
}

void visualixer::sCursorPosition(GLFWwindow * window, double xpos, double ypos){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == window){
			_vInstances.at(i)->onCursorPosition(xpos, ypos);
			return;
		}
	}
	return;
}


//*************************************** PointCloud Visualixer Class *******************************************
cloud_visualixer::cloud_visualixer(){

	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Cloud Visualixer";
	lock_rotation = false;
	//color_ramp = NULL;
	colorby = NULL;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;

	model_centroid[0] = 0.0;
  model_centroid[1] = 0.0;
  model_centroid[2] = 0.0;
}

cloud_visualixer::~cloud_visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}

	delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
}

void cloud_visualixer::set_test_case(){
	num_vertices = 100;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<10; i++){
		for (unsigned int j=0; j<10; j++){
			vertices[(i*10+j)*num_per_vertex] = GLfloat(i);
			vertices[(i*10+j)*num_per_vertex + 1] = GLfloat(j);
			vertices[(i*10+j)*num_per_vertex + 2] = 10.0;
			if (j%2==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 0.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f;
			}
		}
	}

	num_elements = num_vertices;
	num_per_element = 1;
	elements = new GLuint[num_elements*num_per_element];
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	for (unsigned int i=0; i<num_vertices; i++){
		cout << "x: " << vertices[i*num_per_vertex] << " y: " << vertices[i*num_per_vertex +1] << " z: " << vertices[i*num_per_vertex+2] << endl;
	}

	model_centroid[0] = 4.5;
	model_centroid[1] = 4.5;
	model_centroid[2] = 0.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}

void cloud_visualixer::add_cloud(PointCloud * cloud){
	num_vertices = cloud->pointcount;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
		vertices[i*num_per_vertex] = cloud->x[i];
		vertices[i*num_per_vertex + 1] = cloud->y[i];
		vertices[i*num_per_vertex + 2] = cloud->z[i];

	}

	num_elements = num_vertices;
	num_per_element = 1;
	elements = new GLuint[num_elements];
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

  if (cloud->RGB != NULL){
    for (unsigned int i=0; i<cloud->pointcount; i){
      vertices[i*num_per_vertex + 3] = cloud->RGB[i].R/65,535;
      vertices[i*num_per_vertex + 4] = cloud->RGB[i].G/65,535;
      vertices[i*num_per_vertex + 5] = cloud->RGB[i].B/65,535;
    }
  }
  else {
    //rgb ptcolor;
    for (unsigned int i=0; i<cloud->pointcount; i++){
      //ptcolor = color_ramp.get_ramp_color((cloud->z[i]-cloud->zmin)/(cloud->zmax - cloud->zmin));
      vertices[i*num_per_vertex + 3] = 0.0;
      vertices[i*num_per_vertex + 4] = 0.0;
      vertices[i*num_per_vertex + 5] = 1.0;
    }
  }

	model_centroid[0] = (cloud->xmax + cloud->xmin)/2.0;
	model_centroid[1] = (cloud->ymax + cloud->ymin)/2.0;
	model_centroid[2] = (cloud->zmax + cloud->zmin)/2.0;
	xmax = cloud->xmax;
	ymax = cloud->ymax;
	zmax = cloud->zmax;
	xmin = cloud->xmin;
	ymin = cloud->ymin;
	zmin = cloud->zmin;

	// default color by Z
	colorby = new float[cloud->pointcount];
	for (auto i=0; i<cloud->pointcount; i++){
		colorby[i] = cloud->z[i];
	}
	colorby_max = zmax;
	colorby_min = zmin;

	return;
}

bool cloud_visualixer::MainLoop(){

  while(!glfwWindowShouldClose(window_ptr)){

		if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    glfwSwapBuffers(window_ptr);
    glfwPollEvents();

    // Clear the screen to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

    // Draw points
    glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
    //cout << "looping \r" << flush;
	}

	return 0;
}


//*************************************** Mesh Visualixer Class *******************************************
mesh_visualixer::mesh_visualixer(){

	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Mesh Visualixer";
	lock_rotation = false;
	//color_ramp = NULL;
	colorby = NULL;
	vertices = NULL;
	elements = NULL;
	//normals = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;
	num_line_elements = 0;
	//num_normals = 0;

	model_centroid[0] = 0.0;
  model_centroid[1] = 0.0;
  model_centroid[2] = 0.0;
}

mesh_visualixer::~mesh_visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}

	delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	//if (normals != NULL) delete[] normals;
}


void mesh_visualixer::add_mesh(Mesh * mesh){
	Node * node;
	std::map<unsigned int, unsigned int> key_to_index_map;

	num_vertices = mesh->get_num_nodes();
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
		node = mesh->get_node_ptr(mesh->get_node_key(i));

		vertices[i*num_per_vertex] = node->x;
		vertices[i*num_per_vertex + 1] = node->y;
		vertices[i*num_per_vertex + 2] = node->z;
		if (node->boundary){
			vertices[i*num_per_vertex + 3] = 1.0f;
			vertices[i*num_per_vertex + 4] = 0.0f;
			vertices[i*num_per_vertex + 5] = 0.0f;
		}
		else {
			vertices[i*num_per_vertex + 3] = 1.0f;
			vertices[i*num_per_vertex + 4] = 1.0f;
			vertices[i*num_per_vertex + 5] = 1.0f;
		}

		key_to_index_map[node->key] = i;
	}

	// figure out how many line elements are needed
	num_line_elements = 0;
	for (unsigned int i=0; i<num_vertices; i++){
		node = mesh->get_node_ptr(mesh->get_node_key(i));
		num_line_elements+=node->neighbor_keys.size();
	}


	num_elements = num_vertices;
	num_per_element = 1;
	num_per_line_element = 2;
	elements = new GLuint[num_elements*num_per_element + num_line_elements*num_per_line_element];
	line_element_offset = num_elements*num_per_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}

	// DYLAN_TODO: fix this... the key value isn't necessarily the same as the index
	unsigned int elements_added = 0;
	for (unsigned int i=0; i<num_vertices; i++){
		node = mesh->get_node_ptr(mesh->get_node_key(i));

		for (unsigned int j=0; j<node->neighbor_keys.size(); j++){
			//cout << "node: " << node->key << " neighbor: " << node->neighbor_keys[j] << endl;

			elements[line_element_offset + elements_added*num_per_line_element] = key_to_index_map.at(node->key);
			elements[line_element_offset + elements_added*num_per_line_element + 1] = key_to_index_map.at(node->neighbor_keys[j]);

			elements_added++;
		}
	}

	xmax = mesh->get_xmax();
	ymax = mesh->get_ymax();
	zmax = mesh->get_zmax();
	xmin = mesh->get_xmin();
	ymin = mesh->get_ymin();
	zmin = mesh->get_zmin();

	model_centroid[0] = (xmax + xmin)/2.0;
	model_centroid[1] = (ymax + ymin)/2.0;
	model_centroid[2] = (zmax + zmin)/2.0;

	return;
}

void mesh_visualixer::set_test_case(){
	num_vertices = 100;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<10; i++){
		for (unsigned int j=0; j<10; j++){
			vertices[(i*10+j)*num_per_vertex] = GLfloat(i);
			vertices[(i*10+j)*num_per_vertex + 1] = GLfloat(j);
			vertices[(i*10+j)*num_per_vertex + 2] = 10.0;
			if (j%2==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 0.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f;
			}
			else if (j%3==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 0.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f;
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f;
			}
		}
	}

	num_elements = num_vertices;
	num_per_element = 1;
	num_line_elements = 18;
	num_per_line_element = 2;
	elements = new GLuint[num_elements*num_per_element + num_line_elements*num_per_line_element];
	line_element_offset = num_elements*num_per_element;
	// set the point elements
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
	}
	// add in the line elements
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + i*num_per_line_element] = i;
		elements[line_element_offset + i*num_per_line_element+1] = (i+1);
	}
	for (unsigned int i=0; i<9; i++){
		elements[line_element_offset + (i+9)*num_per_line_element] = i*10;
		elements[line_element_offset + (i+9)*num_per_line_element + 1] = (i+1)*10;
	}


	model_centroid[0] = 4.5;
	model_centroid[1] = 4.5;
	model_centroid[2] = 10.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}

void mesh_visualixer::onRender(){
	// Create Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

	// create VBO and copy data to it
  glGenBuffers (1, &vbo);

  // visualizer specific data definitions
  // this one happens to be XYRGB
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

	// Create an element array if necessary
	if (num_elements > 0){
	    glGenBuffers(1, &ebo);
	    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (num_elements * num_per_element + num_line_elements * num_per_line_element) * sizeof(GLuint), elements, GL_STATIC_DRAW);
  }


  //cout << "total elements buffered: " << num_elements * num_per_element + num_line_elements * num_per_line_element << endl;
  // enable point size specification
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  return;
}

bool mesh_visualixer::MainLoop(){

  while(!glfwWindowShouldClose(window_ptr)){

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
    glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
    // Draw lines
    glDrawElements(GL_LINES, num_line_elements*num_per_line_element , GL_UNSIGNED_INT, (void *)(line_element_offset * sizeof(GLuint)));

    //cout << "looping \r" << flush;
	}

	return 0;
}

void mesh_visualixer::onExit(){
	glDeleteProgram(shaderProgram);
  glDeleteShader(fragmentShader);
  glDeleteShader(vertexShader);

  if (num_elements > 0) glDeleteBuffers(1, &ebo);
	if (num_line_elements > 0) glDeleteBuffers(1, &lebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	glfwDestroyWindow( window_ptr );
	glfwTerminate();
	return;
}


//*************************************** Mesh Model Visualixer Class *******************************************
mesh_model_visualixer::mesh_model_visualixer(){
  _vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Mesh Model Visualixer";
	lock_rotation = false;
	//color_ramp = NULL;
	colorby = NULL;
	vertices = NULL;
	elements = NULL;
	normals = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;
	num_normals = 0;

	model_centroid[0] = 0.0;
  model_centroid[1] = 0.0;
  model_centroid[2] = 0.0;
}

mesh_model_visualixer::~mesh_model_visualixer(){
  for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}

	delete[] window_name;
	//if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	if (normals != NULL) delete[] normals;
}


void mesh_model_visualixer::add_model(mesh_model * model){
  // declare vars

  // assign basic info
  num_vertices = model->vertex_count;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
  for (unsigned int i=0; i<num_vertices; i++){
      vertices[i*num_per_vertex] = model->vertices[i*3];
      vertices[i*num_per_vertex+1] = model->vertices[i*3+1];
      vertices[i*num_per_vertex+2]= model->vertices[i*3+2];
      if (i%2==0){
          vertices[i*num_per_vertex + 3] = 0.81f;
          vertices[i*num_per_vertex + 4]  = 0.81f;
          vertices[i*num_per_vertex + 5] = 0.81f;
        }
      else {
          vertices[i*num_per_vertex + 3] = 0.41f;
          vertices[i*num_per_vertex + 4]  = 0.41f;
          vertices[i*num_per_vertex + 5] = 0.41f;
        }
  }

  num_elements = num_vertices/3;
	num_per_element = 3;
  elements = new GLuint[num_elements*num_per_element];
	// set the triangle elements
	for (unsigned int i=0; i<num_elements; i++){
		elements[i*num_per_element] = i*num_per_element;
    elements[i*num_per_element+1] = i*num_per_element+1;
    elements[i*num_per_element+2] = i*num_per_element+2;
	}


  model_centroid[0] = 0.0;
	model_centroid[1] = 0.0;
	model_centroid[2] = 0.0;
	xmax = vertices[0];
	ymax = vertices[1];
	zmax = vertices[2];
	xmin = vertices[0];
	ymin = vertices[1];
	zmin = vertices[2];
  for (unsigned int i=0; i<num_vertices; i++){
    if (vertices[i*num_per_vertex] > xmax) xmax = vertices[i*num_per_vertex];
    if (vertices[i*num_per_vertex] < xmin) xmin = vertices[i*num_per_vertex];
    if (vertices[i*num_per_vertex+1] > ymax) ymax = vertices[i*num_per_vertex+1];
    if (vertices[i*num_per_vertex+1] < ymin) ymin = vertices[i*num_per_vertex+1];
    if (vertices[i*num_per_vertex+2] > zmax) zmax = vertices[i*num_per_vertex+2];
    if (vertices[i*num_per_vertex+2] < zmin) zmin = vertices[i*num_per_vertex+2];

    model_centroid[0] += vertices[i*num_per_vertex];
    model_centroid[1] += vertices[i*num_per_vertex+1];
    model_centroid[2] += vertices[i*num_per_vertex+2];
  }
  model_centroid[0] /= num_vertices;
  model_centroid[1] /= num_vertices;
  model_centroid[2] /= num_vertices;

  return;
}


void mesh_model_visualixer::set_test_case(){
  num_vertices = 100;
	num_per_vertex = 6;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	for (unsigned int i=0; i<num_vertices; i++){
			vertices[i*num_per_vertex] = 0.5*GLfloat(i);
			if (i%2==0) vertices[i*num_per_vertex + 1] = 0.5;
      else vertices[i*num_per_vertex +1] = 0.0;
			if (i%3==1) vertices[i*num_per_vertex + 2] = 0.5;
      else vertices[i*num_per_vertex + 2] = 0.0;
			if (i%2==0){
				vertices[i*num_per_vertex + 3] = 1.0f;
				vertices[i*num_per_vertex + 4]  = 0.0f;
				vertices[i*num_per_vertex + 5] = 0.0f;
			}
			else if (i%3==0){
				vertices[i*num_per_vertex + 3] = 0.0f;
				vertices[i*num_per_vertex + 4]  = 1.0f;
				vertices[i*num_per_vertex + 5] = 0.0f;
			}
			else {
				vertices[i*num_per_vertex + 3] = 0.0f;
				vertices[i*num_per_vertex + 4]  = 0.0f;
				vertices[i*num_per_vertex + 5] = 1.0f;
			}
	}


	num_elements = num_vertices-2;
	num_per_element = 3;
  elements = new GLuint[num_elements*num_per_element];
	// set the triangle elements
	for (unsigned int i=0; i<num_vertices-2; i++){
		elements[i*num_per_element] = i;
    elements[i*num_per_element + 1] = i+1;
    elements[i*num_per_element + 2] = i+2;
	}


	model_centroid[0] = num_vertices*0.5*0.5;
	model_centroid[1] = 0.5*0.5;
	model_centroid[2] = 0.0;
	xmax = num_vertices*0.5;
	ymax = 0.5;
	zmax = 0.0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}

const GLchar * mesh_model_visualixer::VertexShaderSource(){
	// vertex shader sources
	const GLchar* vertexSource =
	    "#version 140\n"
	    "in vec3 position;"
	    "in vec3 color;"
	    "in vec3 normal;"
	    "out vec3 Color;"
	    "uniform mat4 model;"
	    "uniform mat4 view;"
    	"uniform mat4 proj;"
	    "void main() {"
	    "   Color = color;"
      "   gl_PointSize = 2.0;"
	    "   gl_Position = proj*view*model*vec4(position, 1.0);"
	    "}";
	return vertexSource;

}

const GLchar * mesh_model_visualixer::FragmentShaderSource(){
  // fragment shader source
	const GLchar* fragmentSource =
    "#version 140\n"
    "in vec3 Color;"
    "out vec4 outColor;"
    "void main() {"
    "   outColor = vec4(Color, 1.0);"
    "}";
  return fragmentSource;
}

void mesh_model_visualixer::onRender(){
  // Create Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  // create VBO and copy data to it
  glGenBuffers (1, &vbo);

  // visualizer specific data definitions
  // this one happens to be XYRGB
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);

	// Create an element array if necessary
	if (num_elements > 0){
    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, (num_elements * num_per_element + num_line_elements * num_per_line_element) * sizeof(GLuint), elements, GL_STATIC_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (num_elements * num_per_element) * sizeof(GLuint), elements, GL_STATIC_DRAW);
  }


	// create a normal buffer
  	if (num_normals > 0){
	glGenBuffers(1, &normalbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
	glBufferData(GL_ARRAY_BUFFER, num_normals*sizeof(GLfloat)*3,  normals, GL_STATIC_DRAW);
	}
  //cout << "total elements buffered: " << num_elements * num_per_element<< endl;
  // enable point size specification
  //glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  // enable face culling
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  // enable depth test culling
  glEnable(GL_DEPTH_TEST);

  return;
}

void mesh_model_visualixer::onShaders(){

  // create and compile the vertex shader
  vertexShader = glCreateShader(GL_VERTEX_SHADER);
  const GLchar * vssource = this->VertexShaderSource();
  glShaderSource(vertexShader, 1, &(vssource), NULL);
  glCompileShader(vertexShader);

  GLint status;
  GLint logLength;
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE){
		cout << "vertex shader failed to compile" << endl;
		glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH , &logLength);
		if (logLength > 1) {
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(vertexShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}

  // create and compile the fragment shader
  fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  const GLchar * fssource = this->FragmentShaderSource();
  glShaderSource(fragmentShader, 1, &(fssource), NULL);
  glCompileShader(fragmentShader);

  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE){
		cout << "fragment shader failed to compile" << endl;
		glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH , &logLength);
		if (logLength > 1) {
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(fragmentShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}


  // link the vertex and fragment shader into a shader program
  shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glBindFragDataLocation(shaderProgram, 0, "outColor");
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  // Specify the layout of the vertex data
  //cout << "num_vertices: " << num_vertices << " num_vertex_points: " << num_vertex_points << " num_per_vertex: " << num_per_vertex << endl;
  GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
  glEnableVertexAttribArray(posAttrib);
  glVertexAttribPointer(posAttrib, num_vertex_points, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat), 0);

  GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
  glEnableVertexAttribArray(colAttrib);
  glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat), (void*)(num_vertex_points * sizeof(GLfloat)));

  GLint normAttrib = glGetAttribLocation(shaderProgram, "normal");
  glEnableVertexAttribArray(normAttrib);
//glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
	glVertexAttribPointer(
	    normAttrib, 3, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat),                                // stride
	    (void*)0                          // array buffer offset
	);

  rotdeg = 0;
	model = glm::rotate(model, rotdeg, glm::vec3(0.0f, 0.0f, 1.0f)); // angle in radians to suppress some output
	uniModel = glGetUniformLocation(shaderProgram, "model");
	glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

	// Set up view matrix
	zoom_level = 1.0f;
	zoom_scale = 1.0f;
	up_vec = glm::vec3(0.0f, 1.0f, 0.0f);
	// calculate the eye z position so that it can view the whole scene
	if (xmax - xmin > ymax - ymin) eyez_init = (xmax-xmin) + zmax;
	else eyez_init = (ymax - ymin) + zmax;
	//cout << "centroid is: " << model_centroid[0] << ", " << model_centroid[1] << ", " << model_centroid[2] << endl;
	focus_vec = glm::vec3(model_centroid[0], model_centroid[1], model_centroid[2]);
	eye_vec = glm::vec3(model_centroid[0], model_centroid[1], eyez_init);
  view = glm::lookAt(
    eye_vec, // camera position
    focus_vec, // the position to be looking at
    up_vec  // the up vector
  );
  recalcCamera();

  uniView = glGetUniformLocation(shaderProgram, "view");
  glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

  // set up projection matrix
  proj = glm::perspective(0.785f, float(DEFAULT_WIDTH)/float(DEFAULT_HEIGHT), 0.05f, 100000.0f);
  uniProj = glGetUniformLocation(shaderProgram, "proj");
  glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

}

bool mesh_model_visualixer::MainLoop(){

  while(!glfwWindowShouldClose(window_ptr)){

		if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    glfwSwapBuffers(window_ptr);
		glfwPollEvents();

    // Clear the screen to black
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

    // Draw triangles
    //glCullFace(GL_BACK);
    glDrawElements(GL_TRIANGLES, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);

    //cout << "looping \r" << flush;
	}

	return 0;
}

void mesh_model_visualixer::onExit(){
  glDeleteProgram(shaderProgram);
  glDeleteShader(fragmentShader);
  glDeleteShader(vertexShader);

  if (num_elements > 0) glDeleteBuffers(1, &ebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	glfwDestroyWindow( window_ptr );
	glfwTerminate();
	return;
}


//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************

// use cmake to compile
#include "model2mesh.hpp"

int main(int argc, char * argv[]){
	// declare vars

	// test the base class
	visualixer * mywindow = new visualixer();
	//cout << "about to run" << endl;
	mywindow->set_test_case();
	mywindow->run();
	//cout << "about to delete mywindow" << endl;
	delete mywindow;

	//cout << "numero dos" << endl;

	// test the point cloud viewer
	cloud_visualixer * mycvis = new cloud_visualixer();
	//mycvis->set_test_case();
	PointCloud * cloud = PointCloud::read_LAS("../testfiles/ComplexSRSInfo.las");
	//PointCloud * cloud = PointCloud::read_LAS("../testfiles/xyzrgb_manuscript.las");
	//PointCloud * cloud = PointCloud::read_LAS("../testfiles/LAS12_Sample_withRGB_Quick_Terrain_Modeler.las");
	mycvis->add_cloud(cloud);
    mycvis->set_color_ramp(CRamp::DIVERGENT_1);
	mycvis->run();
	delete mycvis;

	// test the mesh viewer
	mesh_visualixer * mymvis = new mesh_visualixer();
	Mesh * mesh = Mesh::create_regular_grid(0.1, (unsigned int)50, (unsigned int)50);//, (unsigned int)30);
	//mymvis->set_test_case();
	mymvis->add_mesh(mesh);
	mymvis->run();
	delete mymvis;

	// test a mesh viewer made from a parametric model
	mesh_visualixer * paravis = new mesh_visualixer();
	parametric_model_2d my_param2;
	Mesh * paramesh;
	my_param2.set_model_name("Display Test Circle");
	my_param2.add_physical_property("Epsilon_rel");
	my_param2.add_physical_property("Mu_rel");
	my_param2.add_material("Air", {1.0, 1.0});
	my_param2.add_material("Dielectric", {5.0, 2.0});
	my_param2.add_object(circle(0.75, vertex_2d(0.0, 0.0), my_param2.get_material("Dielectric")));
	my_param2.print_summary();
	paramesh = build_simple_mesh_2d(&my_param2, 0.1, -1.0, 1.0, -1.0, 1.0, my_param2.get_material("Air"));
	cout << "built the simple mesh" << endl;
	paravis->add_mesh(paramesh);
	paravis->set_color_ramp(CRamp::DIVERGENT_2);
	paravis->set_colorby(paramesh->get_phys_property_ptr("Epsilon_rel"));
	paravis->run();
	delete paravis;
	//delete my_param2;

  // test the mesh model viewer
  mesh_model_visualixer * mymmodvis = new mesh_model_visualixer();
  //mymmodvis->set_test_case();
  mesh_model * mesh_test_model = mesh_model::read_STL("../testfiles/brain-gear.stl");
  mymmodvis->add_model(mesh_test_model);
  mymmodvis->run();
  delete mymmodvis;

	return 0;

}
