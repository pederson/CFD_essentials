#include "./visualixer.hpp"

using namespace std;

// this contains definitions for the openGL visualizer widget
// 
// the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds

vector<visualixer*> _vInstances;

visualixer::visualixer(){

	// really should keep track of instances of visualixer
	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Visualixer";
	lock_rotation = false;
	color_ramp = NULL;
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
	if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
}

char * visualixer::get_window_name(){ return window_name;}

void visualixer::set_window_name(char * w_name){
	window_name = w_name;
	return;
}

void visualixer::run(){
	onInit();
	onRender();
	onShaders();
	MainLoop();
	onExit();
	return;
}

void visualixer::set_test_case(){
	num_vertices = 4;
	num_per_vertex = 5;
	num_vertex_points = 2;
	vertices = new GLfloat[num_vertices*num_per_vertex];
	vertices[0] = 10.0-0.5; vertices[1] = 10.0+0.5; vertices[2] = 1.0; vertices[3] = 0.0; vertices[4] = 0.0;
	vertices[5] = 10.0+0.5; vertices[6] = 10.0+0.5; vertices[7] = 0.0; vertices[8] = 1.0; vertices[9] = 0.0;
	vertices[10] = 10.0+0.5; vertices[11] = 10.0-0.5; vertices[12] = 0.0; vertices[13] = 0.0; vertices[14] = 1.0;
	vertices[15] = 10.0-0.5; vertices[16] = 10.0-0.5; vertices[17] = 1.0; vertices[18] = 1.0; vertices[19] = 1.0;

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
			//cout << "left mouse click " << (state? "on":"off") << endl;
			if (state) {
				glfwGetCursorPos(window_ptr, &x_upon_click, &y_upon_click);
				new_eye = eye_vec;
			}
			else {
				eye_vec = new_eye;
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
			}
			//cout << "right mouse click " << (state? "on":"off")  << endl;
			break;
		otherwise :
			cout << "unknown mouse button" << endl;
	}
	return;
}

void visualixer::onMouseWheel(double xoffset, double yoffset){
	//cout << "MOUSE WHEEL" << endl;
	//cout << "xoffset: " << xoffset << " yoffset: " << yoffset << endl;
	float zoom_scale_new;
	if (yoffset < 0){
		zoom_level -= 0.05;

	}
	else if (yoffset > 0){
		zoom_level += 0.05;
		
	}
	zoom_scale_new = exp(zoom_level);
	eye_vec = focus_vec + (eye_vec-focus_vec)*(zoom_scale_new-zoom_scale + 1.0f);
	zoom_scale = zoom_scale_new;
	view = glm::lookAt(
        eye_vec,
        focus_vec,
        up_vec
    );
}

void visualixer::onKeyboard(int key, int scancode, int action, int modifiers){
	cout << "KEYBOARD PRESS" << endl;
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS){ // resets the view
		eye_vec.x = model_centroid[0];
		eye_vec.y = model_centroid[1];
		eye_vec.z = zmax+1.0f;
		focus_vec.x = model_centroid[0];
		focus_vec.y = model_centroid[1];
		focus_vec.z = model_centroid[2];
		up_vec.x = 0.0f;
		up_vec.y = 1.0f; 
		up_vec.z = 0.0f;

		cout << "zmax is " << zmax << endl;
		view = glm::lookAt(
	        eye_vec,
	        focus_vec,
	        up_vec
	    	);
	}
}

void visualixer::onCursorPosition(double xpos, double ypos){
	//cout << "CURSOR MOVED IN CONTEXT\r" << flush;
	if (left_mouse_engaged) onMouseLeftDrag(xpos, ypos);
	else if(right_mouse_engaged) onMouseRightDrag(xpos, ypos);
	else if(middle_mouse_engaged) onMouseMiddleDrag(xpos, ypos);
}


//DYLAN_TODO: do this using quaternions instead
void visualixer::onMouseLeftDrag(double xpos, double ypos){
	double new_x_pos, new_y_pos;
	int width, height;
	glm::vec3 in_plane, rot_vec;
	glm::vec4 world_in_plane;

	glfwGetWindowSize(window_ptr, &width, &height);
	glfwGetCursorPos(window_ptr, &new_x_pos, &new_y_pos);

	rotdeg = ((new_x_pos-x_upon_click)*(new_x_pos-x_upon_click)
					 + (new_y_pos-y_upon_click)*(new_y_pos-y_upon_click))/(width*width + height*height)*3.14159;
	in_plane = glm::vec3((new_x_pos-x_upon_click), -(new_y_pos-y_upon_click), 0.0f);
	world_in_plane = glm::inverse(proj*view*model)*glm::vec4(in_plane, 1.0f);
	world_in_plane.x = world_in_plane.x - eye_vec.x;
	world_in_plane.y = world_in_plane.y - eye_vec.y;
	world_in_plane.z = world_in_plane.z - eye_vec.z;

	// comes from the cross product
	rot_vec.x = (eye_vec.y-focus_vec.y)*world_in_plane.z - (eye_vec.z-focus_vec.z)*world_in_plane.y;
	rot_vec.y = (eye_vec.z-focus_vec.z)*world_in_plane.x - (eye_vec.x-focus_vec.x)*world_in_plane.z;
	rot_vec.z = (eye_vec.x-focus_vec.x)*world_in_plane.y - (eye_vec.y-focus_vec.y)*world_in_plane.x;


	new_eye = focus_vec +  glm::rotate((eye_vec-focus_vec), rotdeg, -rot_vec);
	view = glm::lookAt(new_eye, focus_vec, up_vec);
}

void visualixer::onMouseRightDrag(double xpos, double ypos){
	//cout << "dragging right mouse button" << endl;
	double new_x_pos, new_y_pos;
	int width, height;
	glm::vec3 in_plane, rot_vec, wip3;
	glm::vec4 world_in_plane;
	GLfloat pan_scale, in_plane_norm;

	glfwGetWindowSize(window_ptr, &width, &height);
	glfwGetCursorPos(window_ptr, &new_x_pos, &new_y_pos);

	
	pan_scale = ((new_x_pos-x_upon_click)*(new_x_pos-x_upon_click) + (new_y_pos-y_upon_click)*(new_y_pos-y_upon_click))/(width*width + height*height)*((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin));
	cout << "pan scale: " << pan_scale << endl;

	in_plane_norm = (new_x_pos-x_upon_click)*(new_x_pos-x_upon_click) + (new_y_pos-y_upon_click)*(new_y_pos-y_upon_click);
	in_plane = glm::vec3((new_x_pos-x_upon_click)/glm::sqrt(in_plane_norm), -(new_y_pos-y_upon_click)/sqrt(in_plane_norm), 0.0f);
	world_in_plane = glm::inverse(proj*view*model)*glm::vec4(in_plane, 1.0f);


	wip3.x = world_in_plane.x;
	wip3.y = world_in_plane.y;
	wip3.z = world_in_plane.z;

	cout << "about to translate " << endl;
	new_eye = eye_vec + pan_scale*wip3;
	new_focus = focus_vec + pan_scale*wip3;
	view = glm::lookAt(eye_vec, focus_vec, up_vec);
	cout << "finished translating" << endl;
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





const GLchar * visualixer::VertexShaderSource(){
	// Shader sources
	cout << "this is the base vertex shader source!" << endl;
	const GLchar* vertexSource =
	    "#version 140\n"
	    "in vec2 position;"
	    "in vec3 color;"
	    "out vec3 Color;"
	    "uniform mat4 model;"
	    "uniform mat4 view;"
    	"uniform mat4 proj;"
	    "void main() {"
	    "   Color = color;"
	    "   gl_Position = proj*view*model*vec4(position, 0.0, 1.0);"
	    "}";
	return vertexSource;

}

const GLchar * visualixer::FragmentShaderSource(){
	cout << "this is the base fragment shader source!" << endl;
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

void visualixer::onRender(){
	// Create Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

	// create VBO and copy data to it
    glGenBuffers (1, &vbo);

    // visualizer specific data definitions
    // this one happens to be XYRGB
    //cout << "binding vbo" << endl;
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);
    
	//cout << "binding ebo" << endl;
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
	//cout << "compiling vertex shader" << endl;

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
		if (logLength > 1)
		{
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(vertexShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}

	//cout << "compiling fragment shader" << endl;

    // create and compile the fragment shader
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    const GLchar * fssource = this->FragmentShaderSource();
    glShaderSource(fragmentShader, 1, &(fssource), NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
	if (status != GL_TRUE){
		cout << "fragment shader failed to compile" << endl;
		glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH , &logLength);
		if (logLength > 1)
		{
		    GLchar* compiler_log = (GLchar*)malloc(logLength);
		    glGetShaderInfoLog(fragmentShader, logLength, 0, compiler_log);
		    printf("%s\n", compiler_log);
		    free (compiler_log);
		}
	}
	

	//cout << "linking shader programs" << endl;
    // link the vertex and fragment shader into a shader program
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glBindFragDataLocation(shaderProgram, 0, "outColor");
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    // Specify the layout of the vertex data
    cout << "num_vertices: " << num_vertices << " num_vertex_points: " << num_vertex_points << " num_per_vertex: " << num_per_vertex << endl;
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

	zoom_level = 0.0f;
	zoom_scale = 1.0f;
	up_vec = glm::vec3(0.0f, 1.0f, 0.0f);
	//cout << "centroid is: " << model_centroid[0] << ", " << model_centroid[1] << ", " << model_centroid[2] << endl;
	focus_vec = glm::vec3(model_centroid[0], model_centroid[1], model_centroid[2]);
	//eye_vec = glm::vec3(model_centroid[0], model_centroid[1], 5.0*(model_centroid[3]+1.0f));
	eye_vec = glm::vec3(model_centroid[0], model_centroid[1], zmax + 1.0f);
    view = glm::lookAt(
        eye_vec, // camera position
        focus_vec, // the position to be looking at
        up_vec  // the up vector
    );
    uniView = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

    // set up projection matrix
    proj = glm::perspective(0.785f, float(DEFAULT_WIDTH)/float(DEFAULT_HEIGHT), 0.05f, 100000.0f);
    //proj = glm::perspective(3.14f/2, 800.0f / 600.0f, 1.0f, 10.0f);
    uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

}


bool visualixer::MainLoop(){

    //cout << "entering main loop" << endl;
    while(!glfwWindowShouldClose(window_ptr)){

		if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    	glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    	//cout << " got keys" << endl;
    	glfwSwapBuffers(window_ptr);
    	//cout << "swapped buffers" << endl;
		glfwPollEvents();
		//cout << "polled events" << endl;

        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        //cout << "cleared colors" << endl;

        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

        // Draw a triangle from the 3 vertices
        glDrawElements(GL_TRIANGLES, num_per_element*num_elements, GL_UNSIGNED_INT, 0);
        //cout << "drew elements" << endl;
        //cout << "looping \r" << flush;
	
		

	}
	//cout << "finished main loop I guess" << endl;

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


/************************************************************/
cloud_visualixer::cloud_visualixer(){

	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Cloud Visualixer";
	lock_rotation = false;
	color_ramp = NULL;
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
	if (color_ramp != NULL) delete[] color_ramp;
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
			vertices[(i*10+j)*num_per_vertex + 2] = 10.0;//GLfloat(i);
			if (j%2==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 4] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
			}
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
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

	model_centroid[0] = 4.5;//num_vertices/2.0;
	model_centroid[1] = 4.5;//num_vertices/2.0;
	model_centroid[2] = 0.0; //num_vertices/2.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0; //num_vertices-1;
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
		if (cloud->RGB != NULL){
			vertices[i*num_per_vertex + 3] = cloud->RGB[i].R/65,535;
			vertices[i*num_per_vertex + 4] = cloud->RGB[i].G/65,535;
			vertices[i*num_per_vertex + 5] = cloud->RGB[i].B/65,535;
		}
		else {
			vertices[i*num_per_vertex + 3] = 1.0;
			vertices[i*num_per_vertex + 4] = 1.0;
			vertices[i*num_per_vertex + 5] = 1.0;
		}
	}

	num_elements = num_vertices;
	num_per_element = 1;
	elements = new GLuint[num_elements];
	for (unsigned int i=0; i<num_vertices; i++){
		elements[i] = i;
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
	return;
}

const GLchar * cloud_visualixer::VertexShaderSource(){
	// Shader sources
	cout << "this is the derived vertex shader source" << endl;
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

bool cloud_visualixer::MainLoop(){

    //cout << "entering main loop" << endl;
    while(!glfwWindowShouldClose(window_ptr)){

		if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    	glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    	//cout << " got keys" << endl;
    	glfwSwapBuffers(window_ptr);
    	//cout << "swapped buffers" << endl;
		glfwPollEvents();
		//cout << "polled events" << endl;

        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        //cout << "cleared colors" << endl;

        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

        // Draw a triangle from the 3 vertices
        glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
        //cout << "drew " << num_elements << " elements" << endl;
        //cout << "looping \r" << flush;
	}
	//cout << "finished main loop I guess" << endl;

	return 0;
}

//*********************************************************************

mesh_visualixer::mesh_visualixer(){

	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Mesh Visualixer";
	lock_rotation = false;
	color_ramp = NULL;
	vertices = NULL;
	elements = NULL;

	num_vertices = 0;
	num_per_vertex = 0;
	num_elements = 0;
	num_line_elements = 0;

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
	if (color_ramp != NULL) delete[] color_ramp;
	if (vertices != NULL) delete[] vertices;
	if (elements != NULL) delete[] elements;
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
		vertices[i*num_per_vertex + 3] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
		vertices[i*num_per_vertex + 4] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
		vertices[i*num_per_vertex + 5] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
		
		key_to_index_map[node->key] = i;
		if (node->key == 1500) cout << "REACHED THE 1500" << endl;
		if (node->key == 1499) cout << "REACHED THE 1499" << endl;
	}
	cout << "yeeerp" << endl;

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
	// set the line elements NOTE: THIS IS WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// DYLAN_TODO: fix this... the key value isn't necessarily the same as the index
	unsigned int elements_added = 0;
	for (unsigned int i=0; i<num_vertices; i++){
		node = mesh->get_node_ptr(mesh->get_node_key(i));

		for (unsigned int j=0; j<node->neighbor_keys.size(); j++){
			//cout << "node: " << node->key << " neighbor: " << node->neighbor_keys[j] << endl;
			//if (node->neighbor_keys[j]==1500) continue;

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
			vertices[(i*10+j)*num_per_vertex + 2] = 10.0;//GLfloat(i);
			if (j%2==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 4] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
			}
			else if (j%3==0){
				vertices[(i*10+j)*num_per_vertex + 3] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 5] = 0.0f; //GLfloat(i)/GLfloat(num_vertices);
			}
			else {
				vertices[(i*10+j)*num_per_vertex + 3] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 4] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
				vertices[(i*10+j)*num_per_vertex + 5] = 1.0f; //GLfloat(i)/GLfloat(num_vertices);
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

	/*
	for (unsigned int i=0; i<18; i++){
		cout << "line element : " << elements[line_element_offset + i*num_per_line_element] << " line element2: " << elements[line_element_offset + i*num_per_line_element+1] << endl;
	}
	*/
	
	//for (unsigned int i=0; i<num_vertices; i++){
	//	cout << "x: " << vertices[i*num_per_vertex] << " y: " << vertices[i*num_per_vertex +1] << " z: " << vertices[i*num_per_vertex+2] << endl;
	//}

	model_centroid[0] = 4.5;//num_vertices/2.0;
	model_centroid[1] = 4.5;//num_vertices/2.0;
	model_centroid[2] = 10.0; //num_vertices/2.0;
	xmax = num_vertices-1;
	ymax = num_vertices-1;
	zmax = 10.0; //num_vertices-1;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	return;
}
const GLchar * mesh_visualixer::VertexShaderSource(){
	// Shader sources
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
	    "   gl_PointSize = 3.0;"
	    "   gl_Position = proj*view*model*vec4(position, 1.0);"
	    "}";
	return vertexSource;
}

void mesh_visualixer::onRender(){
	// Create Vertex Array Object
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

	// create VBO and copy data to it
    glGenBuffers (1, &vbo);

    // visualizer specific data definitions
    // this one happens to be XYRGB
    //cout << "binding vbo" << endl;
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, num_vertices * num_per_vertex * sizeof(GLfloat), vertices, GL_STATIC_DRAW);
    
	//cout << "binding ebo" << endl;
	// Create an element array if necessary
	if (num_elements > 0){
	    glGenBuffers(1, &ebo);
	    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (num_elements * num_per_element + num_line_elements * num_per_line_element) * sizeof(GLuint), elements, GL_STATIC_DRAW);
    }

    cout << "total elements buffered: " << num_elements * num_per_element + num_line_elements * num_per_line_element << endl;
    // enable point size specification
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    return;
}

bool mesh_visualixer::MainLoop(){
	//cout << "entering main loop" << endl;
    while(!glfwWindowShouldClose(window_ptr)){

		if (glfwGetKey(window_ptr, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    	glfwSetWindowShouldClose(window_ptr, GL_TRUE);

    	//cout << " got keys" << endl;
    	glfwSwapBuffers(window_ptr);
    	//cout << "swapped buffers" << endl;
		glfwPollEvents();
		//cout << "polled events" << endl;

        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        //cout << "cleared colors" << endl;

        glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));
        glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

        // Draw nodes
        glDrawElements(GL_POINTS, num_elements*num_per_element , GL_UNSIGNED_INT, NULL);
        // Draw lines
        glDrawElements(GL_LINES, num_line_elements*num_per_line_element , GL_UNSIGNED_INT, (void *)(line_element_offset * sizeof(GLuint)));

        //cout << "drew " << num_elements << " elements" << endl;
        //cout << "looping \r" << flush;
	}
	//cout << "finished main loop I guess" << endl;

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

//*********************************************************************

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
	mycvis->run();
	delete mycvis;

	// test the mesh viewer
	mesh_visualixer * mymvis = new mesh_visualixer();
	Mesh * mesh = Mesh::create_regular_grid(0.1, (unsigned int)50, (unsigned int)30);//, (unsigned int)20);
	//mymvis->set_test_case();
	mymvis->add_mesh(mesh);
	mymvis->run();
	delete mymvis;

	return 0;
	
}