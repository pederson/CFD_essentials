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

	int argc = 1;
	char * argv = "visualixer";

	visualixer_active = false;
	window_name = "Visualixer";
	lock_rotation = false;
	color_ramp = NULL;



	//glutInit(&argc, &argv);
}

visualixer::~visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}
}

char * visualixer::get_window_name(){ return window_name;}

void visualixer::set_window_name(char * w_name){
	window_name = w_name;
	return;
}

void visualixer::run(){
	this->run_test_triangle();
	return;
}

void visualixer::onIdle(){
	cout << "idling \r" << flush;
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

void visualixer::onMouseAction(int button, int updown, int x, int y){
	if (button == MOUSE_SCROLL_UP || button == MOUSE_SCROLL_DOWN ){
		onMouseWheel(button, updown, x, y);
	}
	else if (button == MOUSE_LEFT_BUTTON || button == MOUSE_MIDDLE_BUTTON || button == MOUSE_RIGHT_BUTTON){
		onMouseClick(button, updown, x, y);
	}
	else cout << "Mouse button not recognized" << endl;
}

void visualixer::onMouseClick(int button, int updown, int x, int y){
	bool state;
	if (updown == GLUT_DOWN) state = true;
	else if(updown == GLUT_UP) state = false;
	else cout << "unknown mouse state" << endl;

	switch (button){
		case MOUSE_LEFT_BUTTON :
			left_mouse_engaged = state;
			cout << "left mouse click " << (state? "on":"off") << endl;
			break;
		case MOUSE_MIDDLE_BUTTON :
			middle_mouse_engaged = state;
			cout << "middle mouse click " << (state? "on":"off")  << endl;
			break;
		case MOUSE_RIGHT_BUTTON :
			right_mouse_engaged = state;
			cout << "right mouse click " << (state? "on":"off")  << endl;
			break;
		otherwise :
			cout << "unknown mouse button" << endl;
	}
}

void visualixer::onMouseWheel(int wheel_number, int direction, int x, int y){
	cout << "mouse wheel" << endl;
	//glm::mat4 scalematrix = glm::scale(2.0f, 2.0f, 2.0f);
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
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(2, 0);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(640, 480);
	glut_window_number = glutCreateWindow(window_name);

	// define callbacks (static functions required)
	glutDisplayFunc(sDisplay); //check
	glutReshapeFunc(sReshape); //check
	glutMouseFunc(sMouseAction); //check
	glutMotionFunc(sMouseClickDrag); // check
	//glutMouseWheelFunc(sMouseWheel); // unsure if this works... can live without it
	//glutCloseFunc(sClose);
	glutKeyboardFunc(sKeyDown); // check
	glutKeyboardUpFunc(sKeyUp); // check
	glutSpecialFunc(sSpecialKeyDown); // check // sets special keyboard callback (F and arrow keys)
	glutSpecialUpFunc(sSpecialKeyUp); // check // sets special keyboard release callback
	glutIdleFunc(sIdle); //check
	return;
}

bool visualixer::MainLoop(){
	visualixer_active = true;
	glutMainLoop();
	return true;
}

void visualixer::onExit(){
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
	onInit();
	glutDisplayFunc(test_triangle);
	MainLoop();
	onExit();
	return;
}

void visualixer::sClose(){

}

void visualixer::sReshape(int w, int h){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->glut_window_number == current_window){
			_vInstances.at(i)->onReshape(w, h);
			return;
		}
	}
	return;
}

void visualixer::sDisplay(){

}

void visualixer::sMouseAction(int button, int updown, int x, int y){
	int current_window = glutGetWindow();
	
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->glut_window_number == current_window){
			_vInstances.at(i)->onMouseAction(button, updown, x, y);
			return;
		}
	}
	return;
}

void visualixer::sMouseWheel(int wheel_number, int direction, int x, int y){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->glut_window_number == current_window){
			_vInstances.at(i)->onMouseWheel(wheel_number, direction, x, y);
			return;
		}
	}
	return;
}

void visualixer::sMouseClickDrag(int x, int y){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->glut_window_number == current_window){
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
		if (_vInstances.at(i)->glut_window_number == current_window){
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
		if (_vInstances.at(i)->glut_window_number == current_window){
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
		if (_vInstances.at(i)->glut_window_number == current_window){
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
		if (_vInstances.at(i)->glut_window_number == current_window){
			//_vInstances.at(i)->onSpecialKeyDown(key, x, y);
			return;
		}
	}
	return;
}

void visualixer::sIdle(){
	int current_window = glutGetWindow();

	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i)->glut_window_number == current_window){
			_vInstances.at(i)->onIdle();
			return;
		}
	}
	return;
}

/************************************************************/
cloud_visualixer::cloud_visualixer(){

	_vInstances.push_back(this);

	int argc = 1;
	char * argv = "cloud_visualixer";

	visualixer_active = false;
	window_name = "Cloud Visualixer";
	lock_rotation = false;
	color_ramp = NULL;

	glutInit(&argc, &argv);
}

cloud_visualixer::~cloud_visualixer(){
	for (unsigned int i=0; i<_vInstances.size(); i++){
		if (_vInstances.at(i) == this){
			_vInstances.erase(_vInstances.begin() + i);
			return;
		}
	}
}


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

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	cout << "yup" << endl;

	GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL", NULL, NULL); // Windowed
	//GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL", glfwGetPrimaryMonitor(), nullptr); // Fullscreen

	glfwMakeContextCurrent(window);
	cout << "made context current" << endl;
	while(!glfwWindowShouldClose(window)){
		glfwSwapBuffers(window);
		glfwPollEvents();
		//if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    	//glfwSetWindowShouldClose(window, GL_TRUE);
	}

	//std::this_thread::sleep_for(std::chrono::seconds(1));
	glfwTerminate();
	
	return 0;

	// test the base class
	visualixer * mywindow = new visualixer();
	cout << "about to run" << endl;
	mywindow->run();
	cout << "about to delete mywindow" << endl;
	delete mywindow;

	cout << "numero dos" << endl;

	// test the point cloud viewer
	cloud_visualixer * mycvis = new cloud_visualixer();
	mycvis->run();
	delete mycvis;

	// test the mesh viewer

	return 0;
	
}