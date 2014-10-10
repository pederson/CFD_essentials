/************************************************************************************//**
 * \file function_approx.cpp
 * 
 * File filled with necessary code for function approximation
 *
 ***************************************************************************************/

#include "function_approx.hpp"

//#define _TEST_

#define PI 3.14159265359

FuncApprox::FuncApprox(){

}

FuncApprox::~FuncApprox(){

}


/************************************************************************************//**
 * \brief Determine coefficients of a Chebyshev basis function
 * 
 *  This function calculates the coefficients of Chebyshev basis functions
 *  where bk = -cos(k*x)
 *
 *  \param function_handle : handle to real function to be approximated
 *  \param M : the number of basis functions to use
 *
 ***************************************************************************************/
double * FuncApprox::basis_cheby(double (*function_handle)(double), unsigned int M){
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

/************************************************************************************//**
 * \brief Determine coefficients of a piecewise linear basis using lagrange polynomials
 * 
 *  This function calculates the coefficients of linear lagrange polynomical basis functions
 *  where bk = delta_function(xk)
 *
 *  \param function_handle : handle to real function to be approximated
 *  \param M : the number of basis functions to use
 *
 ***************************************************************************************/
double * FuncApprox::basis_lagrange_poly(double (*function_handle)(double), unsigned int M){
	// declare vars
	unsigned int *ak;
	double * xi;

	// initialize all x values
	xi = new double[M+1];
	for (unsigned int i=0; i<M+1; i++) xi[i] = i*2.0/M - 1.0;

	// the ak values are just the values of the function at the standard points
	ak = new double[M+1];
	for (unsigned int i=0; i<M+1; i++) ak[i] = function_handle(xi[i]);

	// delete stuff
	delete[] xi;

	return ak;
}

/************************************************************************************//**
 * \brief Evaluate approximate values using a Chebyshev function approximation
 * 
 *  This function takes x values from [-1, 1] and approximates the function values at 
 *  these points by using Chebyshev basis functions. 
 *
 *  \param function_handle : handle to real function to be approximated
 *  \param M : the number of basis functions to use
 *
 ***************************************************************************************/
double * FuncApprox::cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts){
	// declare vars
	unsigned int *ki;
	double *ak, *values;

	// initialize stuff
	values = new double[N_eval_pts];
	ki = new unsigned int[M+1];
	for (unsigned int i=0; i<M+1; i++) ki[i] = i;

	// get the ak coefficients
	ak = basis_cheby(function_handle, M);
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
