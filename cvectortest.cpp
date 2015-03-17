#include "cvector.hpp"

using namespace std;

// compile with:
// 				g++ -std=c++11 cvectortest.cpp -o cvectortest

int main(int argc, char * argv[]){

	//
	unsigned int nints = 20;

	// sequence 
	double integers[nints];

	for (auto i=0; i<nints; i++) integers[i] = i;

	cvector squares, squaresp1, cubes, quads;

	squares.multiply(integers);
	squares.multiply(integers);

	squaresp1.multiply(integers);
	squaresp1.multiply(integers);
	squaresp1.multiply(0.5);
	squaresp1.set_additive_constant(1);

	cubes.multiply(integers);
	cubes.multiply(integers);
	cubes.multiply(integers);

	quads.multiply(integers);
	quads.multiply(integers);
	quads.multiply(integers);
	quads.multiply(integers);

	for (auto i=0; i<nints; i++){
		cout << "orig: " << integers[i] << "  squares: " << squares[i] << "  squares+1: " << squaresp1[i] << "  cubes: " << cubes[i] << "  quads: " << quads[i] << endl;
	}
}