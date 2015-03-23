#include "cvector.hpp"
#include <time.h>

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

	// test the speed of cvector versus raw multiply/divide
	time_t tstart, tstop;
	clock_t cstart, cstop;
	double ttot;
	unsigned int nlots = 200000;
	double * numbers = new double[nlots];
	double * results = new double[nlots];
	cvector mycvec;
	for (auto i=0; i<nlots; i++) numbers[i] = i/1000;

	mycvec.multiply(numbers);
	mycvec.multiply(numbers);
	mycvec.multiply(numbers);
	mycvec.multiply(10);
	mycvec.divide(numbers);
	mycvec.set_additive_constant(1);
	cstart = clock();
	for (auto i=0; i<nlots; i++){
	results[i] = mycvec.at(i);
	}
	cstop = clock();
	ttot = (double) (cstop-cstart) / CLOCKS_PER_SEC * 1000.0;


	cout << "cvector multiplication took " << ttot << " milliseconds. That's " << nlots/ttot << " operations per second" << endl;


	cstart = clock();
	for (auto i=0; i<nlots; i++){
		results[i] = 1 + numbers[i] * numbers[i] * numbers[i] * 10 / numbers[i];
		//cout << i << endl;
	}
	cstop = clock();
	ttot = (double) (cstop-cstart) / CLOCKS_PER_SEC * 1000.0;

	cout << "raw multiplication took " << ttot << " milliseconds. That's " << nlots/ttot << " operations per second" << endl;

	
	delete[] numbers;
	delete[] results;
}