#include "Visualixer.hpp"

using namespace std;

//#define _TEST_

// this contains definitions for the openGL visualizer widget
//
// the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds

std::vector<visualixer*> _vInstances;

//*************************************** Base Visualixer Class *******************************************
visualixer::visualixer(){

	// really should keep track of instances of visualixer
	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Visualixer";
	rotation_lock = false;
	_colorby = nullptr;
	_color_alpha = nullptr;
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

	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
	if (_color_alpha != nullptr) delete[] _color_alpha;
	if (_colorby != nullptr) delete[] _colorby;
}

void visualixer::set_window_name(string w_name){
	window_name = w_name;
	return;
}

void visualixer::set_color_ramp(CRamp ramp_name){
  color_ramp.set_ramp(ramp_name);
}





void visualixer::run(){
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

void visualixer::set_test_case(){
	num_vertices = 4;
	num_per_vertex = 7;
	num_vertex_points = 3;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	vertices[0] = 10.0-0.5; vertices[1] = 10.0+0.5; vertices[2] = 0.0; vertices[3] = 1.0; vertices[4] = 0.0; vertices[5] = 0.0; vertices[6] = 1.0;
	vertices[7] = 10.0+0.5; vertices[8] = 10.0+0.5; vertices[9] = 0.0; vertices[10] = 0.0; vertices[11] = 1.0; vertices[12] = 0.0; vertices[13] = 1.0;
	vertices[14] = 10.0+0.5; vertices[15] = 10.0-0.5; vertices[16] = 0.0; vertices[17] = 0.0; vertices[18] = 0.0; vertices[19] = 1.0; vertices[20] = 0.1;
	vertices[21] = 10.0-0.5; vertices[22] = 10.0-0.5; vertices[23] = 0.0; vertices[24] = 1.0; vertices[25] = 1.0; vertices[26] = 1.0; vertices[27] = 0.1;

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
	if (zoom_level > 5.0) zoom_level = 5.0; // limit the zoom in
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
	else if (key == GLFW_KEY_C && action == GLFW_PRESS){
		cycleColorRamp();
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
	//cout << "dragging middle mouse button" << endl;
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

void visualixer::cycleColorRamp(){
	if (_colorby==nullptr) return;
	color_ramp.cycle_ramp();
	onColors();
	onRender();
	onShaders();
	return;
}





const GLchar * visualixer::VertexShaderSource(){
	// vertex shader sources
	const GLchar* vertexSource =
	    "#version 140\n"
	    "in vec3 position;"
	    "in vec4 color;"
	    "out vec4 Color;"
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
    "in vec4 Color;"
    "out vec4 outColor;"
    "void main() {"
    "   outColor = Color;"//vec4(Color, 1.0);"
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


	window_ptr = glfwCreateWindow(DEFAULT_WIDTH, DEFAULT_HEIGHT, window_name.c_str(), NULL, NULL); // Windowed
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


void visualixer::onPrepareData(){


}


void visualixer::onColors(){
	if (_colorby == nullptr) return;
	if (_colorby_max - _colorby_min == 0.0) return;
	// modify the vertex array to incorporate user-defined colors


	rgb ptcolor;
	for (auto i=0; i<num_vertices; i++){
		ptcolor = color_ramp.get_ramp_color(float((_colorby[i])/(_colorby_max - _colorby_min)));
		vertices[(i+1)*num_per_vertex - 4] = ptcolor.R;
		vertices[(i+1)*num_per_vertex - 3] = ptcolor.G;
		vertices[(i+1)*num_per_vertex - 2] = ptcolor.B;
	}

	return;
}


void visualixer::onAlpha(){
	if (_color_alpha == nullptr) return;

	
	for (auto i=0; i<num_vertices; i++){
		vertices[(i+1)*num_per_vertex - 1] = _color_alpha[i];
	}

}


void visualixer::onRender(){
	// Create Vertex Array Object
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

	// create VBO and copy data to it
  glGenBuffers (1, &vbo);

  // visualizer specific data definitions
  // this one happens to be XYZRGB
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

  // enable alpha channel
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

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
  glVertexAttribPointer(colAttrib, 4, GL_FLOAT, GL_FALSE, num_per_vertex * sizeof(GLfloat), (void*)(num_vertex_points * sizeof(GLfloat)));

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
  //proj = glm::perspective(0.785f, float(DEFAULT_WIDTH)/float(DEFAULT_HEIGHT), 0.05f, 100000.0f);
  proj = glm::perspective(0.785f, float(DEFAULT_WIDTH)/float(DEFAULT_HEIGHT), 0.000005f, 100000.0f);
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





//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************
//*************************************** TEST SECTION *******************************************

#ifdef _TEST_
// use cmake to compile
#include "Converter.hpp"

int main(int argc, char * argv[]){
	// declare vars

	// test the base class
	visualixer * mywindow = new visualixer();
	//cout << "about to run" << endl;
	mywindow->set_test_case();
	mywindow->run();
	//cout << "about to delete mywindow" << endl;
	delete mywindow;

	return 0;

}

#endif
