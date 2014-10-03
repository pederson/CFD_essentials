#include "function_approx.hpp"

#define _TEST_

#define PI 3.14159265359

double * cheby_basis(double (*function_handle)(double), unsigned int M){
	// declare vars
	unsigned int *ki;
	double * si, *ak;
	double sgn;

	// initialize all s values
	si = new double[M+1];
	for (unsigned int i=0; i<M+1; i++) si[i] = i*PI/M;

	// initialize all k values
	ki = new unsigned int[M];
	for (unsigned int i=0; i<M+1; i++) ki[i] = i;

	for (unsigned int i=0; i<M+1; i++) std::cout << "si: " << si[i] << " f(-c(si)): " << function_handle(-cos(si[i])) << std::endl;
	std::cout << "first term: " << -(function_handle(-1.0) + function_handle(1.0))/2.0 << std::endl;
	// evaluate the ak coefficients by the trapezoid rule
	ak = new double[M+1];
	for (unsigned int i=0; i<M+1; i++) ak[i] = 0.0;
	for (unsigned int i=0; i<M+1; i++){ // loop over ki
		// determine the sign of the q(pi)=f(1)
		if (ki[i]%2 == 0) sgn = 1.0;
		else sgn = -1.0;
		ak[i] = -(function_handle(-1.0) + sgn*function_handle(1.0))/2.0;
		for (unsigned int j=1; j<M+1; j++){ // loop over si
			ak[i] += function_handle(-cos(si[j]))*cos(ki[i]*si[j]);
		}
		ak[i] *= 2.0/double(M);
	}

	// delete stuff
	delete[] si;
	delete[] ki;

	return ak;
}

double * approx_cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts){
	// declare vars
	unsigned int *ki;
	double *ak, *values;

	// initialize stuff
	values = new double[N_eval_pts];
	ki = new unsigned int[M+1];
	for (unsigned int i=0; i<M+1; i++) ki[i] = i;

	// get the ak coefficients
	ak = cheby_basis(function_handle, M);
	for (unsigned int i=0; i<M+1; i++) std::cout << "ak[i]: " << ak[i] << std::endl;
	
	// evaluate input points
	for (unsigned int i=0; i<N_eval_pts; i++){
		values[i] = ak[0]/2.0;
		for (unsigned int j=1; j<M+1; j++){
			values[i] += ak[j]*cos(ki[j]*acos(-eval_points[i]));
		}
	}			

	//delete stuff
	delete[] ak;
	delete[] ki;

	return values;
}

#ifdef _TEST_

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;

double gauss_func(double x);
double gauss_func(double x){
	return exp(-x*x/2.0);
}


int main(int argc, char *argv[]){
	// test on a gaussian function
	unsigned int npts = 100, M;
	double * randx = new double[npts];
	double * approxvals;

	// generate random numbers from -1 to 1 to test on 
	srand(time(NULL));
	for (unsigned int i=0; i<npts; i++) randx[i] = double(rand()%1000)/500.0 - 1.0;

	// find the approximate values
	M = 32;
	approxvals = approx_cheby(gauss_func, M, randx, npts);

	for (unsigned int i=0; i<npts; i++){
		cout << "x: " << randx[i] << " gauss value: " << gauss_func(randx[i]) << " approx value: " << approxvals[i] << endl;
	}

	return 0;
}
 
#endif
