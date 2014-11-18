#include "visualixer.hpp"

using namespace std;

// this contains definitions for the openGL visualizer widget
// basically just a placeholder and reminder for now,
// but the visualizer should be able to:
//  - visualize flow in 2d and 3d and move around in it
//  - visualize the geometry in 2d and 3d and move around in it
//  - visualize and interact with point clouds

visualixer::visualixer(){

}
visualixer::~visualixer(){

}
void visualixer::draw_test_triangle(){
	/*
	if ( !glfwInit() ){
		fprintf(stderr, "Failed to initialize GLFW\n");
		throw -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL
	 
	// Open a window and create its OpenGL context
	GLFWwindow* window; // (In the accompanying source code, this variable is global)
	window = glfwCreateWindow( 1024, 768, "Tutorial 01", NULL, NULL);
	if( window == NULL ){
	    fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
	    glfwTerminate();
	    return -1;
	}
	glfwMakeContextCurrent(window); // Initialize GLEW
	glewExperimental=true; // Needed in core profile
	if (glewInit() != GLEW_OK) {
	    fprintf(stderr, "Failed to initialize GLEW\n");
	    return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	 
	do{
	    // Draw nothing, see you in tutorial 2 !
	 
	    // Swap buffers
	    glfwSwapBuffers(window);
	    glfwPollEvents();
	 
	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
	glfwWindowShouldClose(window) == 0 );
	*/

	/*
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT); 
	glColor3f(1.0, 1.0, 1.0);
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	glBegin(GL_POLYGON);
	glVertex2f(-0.5, -0.5);
	glVertex2f(-0.5, 0.5);
	glVertex2f(0.5, 0.5);
	glVertex2f(0.5, -0.5);
	glEnd();
	glFlush();


	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(500,500);
	glutInitWindowPosition(100,100);
	glutCreateWindow("OpenGL - First window demo");
	glutDisplayFunc(renderFunction);
	glutMainLoop();
	return 0; 
	*/

}

int init_resources(void)
{
	GLint compile_ok = GL_FALSE, link_ok = GL_FALSE;

	// Vertex shader ( sets the points)
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	const char *vs_source = 
	#ifdef GL_ES_VERSION_2_0
	"#version 100\n"  // OpenGL ES 2.0
	#else
	"#version 120\n"  // OpenGL 2.1
	#endif
	"attribute vec2 coord2d;                  "
	"void main(void) {                        "
	"  gl_Position = vec4(coord2d, 0.0, 1.0); "
	"}";
	glShaderSource(vs, 1, &vs_source, NULL);
	glCompileShader(vs);
	glGetShaderiv(vs, GL_COMPILE_STATUS, &compile_ok);
	if (0 == compile_ok)
	{
		fprintf(stderr, "Error in vertex shader\n");
		return 0;
	}

	// Fragment shader (fills in the empty space between points)
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
	const char *fs_source =
	"#version 120           \n"
	"void main(void) {        "
	"  gl_FragColor[0] = 0.0; "
	"  gl_FragColor[1] = 0.0; "
	"  gl_FragColor[2] = 1.0; "
	"}";
	glShaderSource(fs, 1, &fs_source, NULL);
	glCompileShader(fs);
	glGetShaderiv(fs, GL_COMPILE_STATUS, &compile_ok);
	if (!compile_ok) {
		fprintf(stderr, "Error in fragment shader\n");
	return 0;
	}

	// GLSL program
	program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);
	glGetProgramiv(program, GL_LINK_STATUS, &link_ok);
	if (!link_ok) {
		fprintf(stderr, "glLinkProgram:");
	return 0;
	}

	// pass the triangle vertices
	const char* attribute_name = "coord2d";
	attribute_coord2d = glGetAttribLocation(program, attribute_name);
	if (attribute_coord2d == -1) {
		fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
	return 0;
	}
  /* FILLED IN LATER */
  return 1;
}
 
void onDisplay()
{
	/* Clear the background as white */
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glUseProgram(program);
	glEnableVertexAttribArray(attribute_coord2d);
	GLfloat triangle_vertices[] = {
	 0.0,  0.8,
	-0.8, -0.8,
	 0.8, -0.8,
	};
	/* Describe our vertices array to OpenGL (it can't guess its format automatically) */
	glVertexAttribPointer(
	attribute_coord2d, // attribute
	2,                 // number of elements per vertex, here (x,y)
	GL_FLOAT,          // the type of each element
	GL_FALSE,          // take our values as-is
	0,                 // no extra data between each position
	triangle_vertices  // pointer to the C array
	);

	/* Push each element in buffer_vertices to the vertex shader */
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glDisableVertexAttribArray(attribute_coord2d);

	/* Display the result */
	glutSwapBuffers();
  /* FILLED IN LATER */
}
 
void free_resources()
{
	glDeleteProgram(program);
  /* FILLED IN LATER */
}

int main(int argc, char * argv[]){
	// declare vars
	//visualixer * glwindow = new visualixer();

	/*
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
		glutCreateWindow("red 3D lighted cube");
		printf("GL_VERSION = %s\n",glGetString(GL_VERSION) ); // GL_VERSION = 2.1.2 NVIDIA 195.36.24 
		return 0; 
	*/

	/*
	char *version = NULL;
	char *vendor = NULL;
	char *renderer = NULL;
	char *extensions = NULL;
	GLuint idWindow = 0;
	int	glutVersion;

	glutInit(&argc, argv);
	glutInitWindowSize(1,1);
	glutInitDisplayMode(GLUT_RGBA);
	idWindow = glutCreateWindow(PROGRAM);
	glutHideWindow();

	glutVersion = glutGet(0x01FC);
	version =     (char*)glGetString(GL_VERSION);
	vendor =      (char*)glGetString(GL_VENDOR);
	renderer =    (char*)glGetString(GL_RENDERER);
	extensions =  (char*)glGetString(GL_EXTENSIONS);

	printf("GLUT=%d\nVERSION=%s\nVENDOR=%s\nRENDERER=%s\nEXTENSIONS=%s\n",
	glutVersion,version,vendor,renderer,extensions);

	glutDestroyWindow(idWindow);
	return(0);
	*/
	
	
	glutInit(&argc, argv);
	glutInitContextVersion(2,0);
	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);
	glutInitWindowSize(640, 480);

	// Extension wrangler initialising //
	glewExperimental = GL_TRUE;
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
		glutDisplayFunc(onDisplay);
		glutMainLoop();
	}

	// If the program exits in the usual way,
	// free resources and exit with a success 
	free_resources();
	return EXIT_SUCCESS;
	
	
}