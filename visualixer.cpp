#include "visualixer.hpp"

using namespace std;

// this contains definitions for the openGL visualizer widget
// basically just a placeholder and reminder for now,
// but the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds

visualixer::visualixer(){

	// really should keep track of instances of visualixer


	int argc = 1;
	char * argv = "visualixer";

	visualixer_active = false;
	window_name = "Visualixer";
	lock_rotation = false;
	color_ramp = NULL;

	glutInit(&argc, &argv);
}

visualixer::~visualixer(){

}

char * visualixer::get_window_name(){ return window_name;}

void visualixer::set_window_name(char * w_name){
	window_name = w_name;
	return;
}

void visualixer::onInit(){
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(2, 0);
	glutInitWindowPosition(200, 100);
	glutInitWindowSize(640, 480);
	glutCreateWindow(window_name);

	// define callbacks
	glutDisplayFunc(sDisplay);
	glutReshapeFunc(sReshape);
	glutMouseFunc(sMouse);
	glutMotionFunc(sMotion);
	glutMouseWheelFunc(sMouseWheel);
	glutCloseFunc(sClose);
	glutKeyboardFunc(sKeyDown);
	glutKeyboardUpFunc(sKeyUp);
	glutIdleFunc(sIdle);
	return;
}

bool visualixer::MainLoop(){
	visualixer_active = true;
	glutMainLoop();
	return true;
}

void visualixer::onExit(){

}

void visualixer::OnKeyDown(int new_key, char cAscii){
	if (cAscii == 27){ // ESC
		cout << "teehee" << endl;
		//this->Close();
	}
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

void visualixer::sClose(void){

}

void visualixer::sReshape(int w, int h){

}

void visualixer::sDisplay(void){

}

void visualixer::sMouse(int button, int updown, int x, int y){

}

void visualixer::sMouseWheel(int wheel_number, int direction, int x, int y){

}

void visualixer::sMotion(int x, int y){

}

void visualixer::sKeyUp(unsigned char key, int x, int y){

}

void visualixer::sKeyDown(unsigned char key, int x, int y){

}

void visualixer::sIdle(void){

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
	visualixer * mywindow = new visualixer();

	// test the base class
	mywindow->run_test_triangle();

	// test the mesh viewer

	// test the point cloud viewer

	return 0;
	
	/*
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(2, 0);
	glutInitWindowPosition(200, 100);
	glutInitWindowSize(640, 480);
	glutCreateWindow("Dylan's window");
	*/

	/*
	// register callbacks
	glutDisplayFunc(renderScene);
	// enter GLUT event processing cycle
	glutMainLoop();
	return 1;
	*/
	
	/*
	// Extension wrangler initialising //
	//glewExperimental = GL_TRUE;
	GLenum glew_status = glewInit();
	if ( glew_status != GLEW_OK)
	{
		//fprintf(stderr, "There was an error\n" );
		fprintf(stderr, "Error in glewInit: %s\n", glewGetErrorString(glew_status));
		return EXIT_FAILURE;
	}

	// When all init functions run without errors,
	//the program can initialise the resources //
	if (init_resources())
	{
		// We can display it if everything goes OK 
		//glutDisplayFunc(onDisplay);
		cout << "im rendering" << endl;
		glutDisplayFunc(renderScene);
		glutMainLoop();
	}
	cout << "i didn't render" << flush;

	// If the program exits in the usual way,
	// free resources and exit with a success 
	free_resources();
	return EXIT_SUCCESS;
	*/
}