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
}

visualixer::~visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}

	delete[] window_name;
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








void visualixer::onReshape(int new_width, int new_height){
	cout << "reshaping" << endl;
	/*
	if (new_height == 0) new_height = 1;
	float ratio = 1.0*new_width/new_height;

	// projection
	glMatrixMode(GL_PROJECTION);

	// reset matrix
	glLoadIdentity();

	// set the viewport to be the whole window
	glViewport(0, 0, new_width, new_height);

	// set the correct perspective
	gluPerspective(45, ratio, 1, 1000);

	// go back to modelview
	glMatrixMode(GL_MODELVIEW);
	*/
}












void visualixer::onMouseClick(int button, int action, int modifiers){
	bool state;
	if (action == GLFW_PRESS) state = true;
	else if(action == GLFW_RELEASE) state = false;
	else cout << "unknown mouse state" << endl;


	switch (button){
		case GLFW_MOUSE_BUTTON_LEFT :
			left_mouse_engaged = state;
			cout << "left mouse click " << (state? "on":"off") << endl;
			glfwGetCursorPos(window_ptr, &x_upon_click, &y_upon_click);
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE :
			middle_mouse_engaged = state;
			cout << "middle mouse click " << (state? "on":"off")  << endl;
			break;
		case GLFW_MOUSE_BUTTON_RIGHT :
			right_mouse_engaged = state;
			cout << "right mouse click " << (state? "on":"off")  << endl;
			break;
		otherwise :
			cout << "unknown mouse button" << endl;
	}
	return;
}

void visualixer::onMouseWheel(double xoffset, double yoffset){
	cout << "MOUSE WHEEL" << endl;
	camZ += 0.05;
	view = glm::lookAt(
        glm::vec3(0.0f, 0.0f, camZ),
        glm::vec3(0.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f)
    );
	//rotdeg = 5*3.14159/180;
	//model = glm::rotate(model, rotdeg, glm::vec3(0.0f, 1.0f, 1.0f));
	//glm::mat4 scalematrix = glm::scale(2.0f, 2.0f, 2.0f);
}

void visualixer::onKeyboard(int key, int scancode, int action, int modifiers){
	cout << "KEYBOARD PRESS" << endl;
}

void visualixer::onCursorPosition(double xpos, double ypos){
	cout << "CURSOR MOVED IN CONTEXT\r" << flush;
	if (left_mouse_engaged){
		double new_x_pos, new_y_pos;
		glfwGetCursorPos(window_ptr, &new_x_pos, &new_y_pos);
		rotdeg = 1*3.14159/180;
		model = glm::rotate(model, rotdeg, glm::vec3(-(new_x_pos-x_upon_click), -(new_y_pos-y_upon_click), 0.0f));
	}
}











void visualixer::onMouseClickDrag(int x, int y){
	if (left_mouse_engaged) onMouseLeftDrag(x, y);
	else if (right_mouse_engaged) onMouseRightDrag(x, y);
	else cout << "dragging middle mouse button" << endl;
}

void visualixer::onMouseLeftDrag(int x, int y){
	cout << "dragging left mouse button" << endl;
	//glm::vec3 rotation_axis(0, 0, 1);
	//glm::rotate(2, rotation_axis);
}

void visualixer::onMouseRightDrag(int x, int y){
	cout << "dragging right mouse button" << endl;
	//glm::mat4 transmatrix = glm::translate(0.5f, 0.0f, 0.0f);
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


	window_ptr = glfwCreateWindow(400, 400, window_name, NULL, NULL); // Windowed
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
    GLfloat vertices[] = {
        -0.5f,  0.5f, 1.0f, 0.0f, 0.0f, // Top-left
         0.5f,  0.5f, 0.0f, 1.0f, 0.0f, // Top-right
         0.5f, -0.5f, 0.0f, 0.0f, 1.0f, // Bottom-right
        -0.5f, -0.5f, 1.0f, 1.0f, 1.0f  // Bottom-left
    };
    //cout << "binding vbo" << endl;

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    
	//cout << "binding ebo" << endl;
	// Create an element array
    glGenBuffers(1, &ebo);
 
    GLuint elements[] = {
        0, 1, 2,
        2, 3, 0
    };
 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
    return;
}

void visualixer::onShaders(){
	//cout << "compiling vertex shader" << endl;

    // create and compile the vertex shader
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    const GLchar * vssource = VertexShaderSource();
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
    const GLchar * fssource = FragmentShaderSource();
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
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), 0);

    GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
    glEnableVertexAttribArray(colAttrib);
    glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));

    rotdeg = 0;
	model = glm::rotate(model, rotdeg, glm::vec3(0.0f, 0.0f, 1.0f)); // angle in radians to suppress some output
	uniModel = glGetUniformLocation(shaderProgram, "model");
	glUniformMatrix4fv(uniModel, 1, GL_FALSE, glm::value_ptr(model));

	// Set up view matrix
	camZ = 1.05f;
    view = glm::lookAt(
        glm::vec3(0.0f, 0.0f, camZ), // camera position
        glm::vec3(0.0f, 0.0f, 0.0f), // the position to be looking at
        glm::vec3(0.0f, 1.0f, 0.0f)  // the up vector
    );
    uniView = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

    // set up projection matrix
    proj = glm::perspective(0.785f, 800.0f / 600.0f, 1.0f, 10.0f);
    //proj = glm::perspective(3.14f/2, 800.0f / 600.0f, 1.0f, 10.0f);
    uniProj = glGetUniformLocation(shaderProgram, "proj");
    glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

}

const GLchar * visualixer::VertexShaderSource(){
	// Shader sources
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
	const GLchar* fragmentSource =
    "#version 140\n"
    "in vec3 Color;"
    "out vec4 outColor;"
    "void main() {"
    "   outColor = vec4(Color, 1.0);"
    "}";
    return fragmentSource;
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
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
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

    glDeleteBuffers(1, &ebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	glfwDestroyWindow( window_ptr );
	glfwTerminate();
	return;
}


void visualixer::onKeyDown(unsigned char key, int x, int y){
	cout << "Pressed Key with unsigned char : " << key << endl;
	return;
}

void visualixer::onKeyUp(unsigned char key, int x, int y){
	cout << "Released Key with unsigned char: " << key << endl;
	return;
}

void visualixer::SetFullscreen(bool bFullscreen){
	if (bFullscreen){
		glutFullScreen();
	}
	else{
		glutPositionWindow(DEFAULT_CENTER_X, DEFAULT_CENTER_Y);
		glutReshapeWindow(DEFAULT_WIDTH, DEFAULT_HEIGHT);
	}
	return;
}

void visualixer::Hide(){
	glutHideWindow();
	return;
}

void visualixer::Show(){
	glutShowWindow();
	return;
}
void visualixer::Close(){
	glutDestroyWindow(0);
	return;
}

void visualixer::run_test_triangle(){
	return;
}

/*
void visualixer::sClose(){

}

void visualixer::sReshape(int w, int h){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			_vInstances.at(i)->onReshape(w, h);
			return;
		}
	}
	return;
}

void visualixer::sDisplay(){

}
*/

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

/*
void visualixer::sMouseClickDrag(int x, int y){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			_vInstances.at(i)->onMouseClickDrag(x, y);
			return;
		}
	}
	return;
}

void visualixer::sKeyUp(unsigned char key, int x, int y){
	int current_window = glutGetWindow();
	//cout << "Pressed Key with unsigned char: " << key << endl;
	
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			_vInstances.at(i)->onKeyUp(key, x, y);
			return;
		}
	}
	return;
}

void visualixer::sKeyDown(unsigned char key, int x, int y){
	int current_window = glutGetWindow();
	//cout << "Pressed Key with unsigned char: " << key << endl;

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			_vInstances.at(i)->onKeyDown(key, x, y);
			return;
		}
	}
	return;
}

void visualixer::sSpecialKeyUp(int key, int x, int y){
	int current_window = glutGetWindow();
	cout << "Released Special Key " << endl;

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			//_vInstances.at(i)->onSpecialKeyUp(key, x, y);
			return;
		}
	}
	return;
}

void visualixer::sSpecialKeyDown(int key, int x, int y){
	int current_window = glutGetWindow();
	cout << "Pressed Special Key " << endl;

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			//_vInstances.at(i)->onSpecialKeyDown(key, x, y);
			return;
		}
	}
	return;
}

void visualixer::sIdle(){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->window_ptr == current_window){
			_vInstances.at(i)->onIdle();
			return;
		}
	}
	return;
}
*/

/************************************************************/
cloud_visualixer::cloud_visualixer(){

	_vInstances.push_back(this);

	visualixer_active = false;
	window_name = "Cloud Visualixer";
	lock_rotation = false;
	color_ramp = NULL;
}

cloud_visualixer::~cloud_visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}
}

// this is old deprecated style of opengl
void test_triangle(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glBegin(GL_TRIANGLES);
		glVertex3f(-0.5,-0.5,0.0);
		glVertex3f(0.5,0.0,0.0);
		glVertex3f(0.0,0.5,0.0);
	glEnd();

    glutSwapBuffers();
}



int main(int argc, char * argv[]){
	// declare vars

	// test the base class
	visualixer * mywindow = new visualixer();
	//cout << "about to run" << endl;
	mywindow->run();
	//cout << "about to delete mywindow" << endl;
	delete mywindow;

	//cout << "numero dos" << endl;

	// test the point cloud viewer
	cloud_visualixer * mycvis = new cloud_visualixer();
	mycvis->run();
	delete mycvis;

	// test the mesh viewer

	return 0;
	
}