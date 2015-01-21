
#include "integration.hpp"

Integrator::Integrator{

}

double Integrator::trapezoid(double (*function_handle)(double), double x_start, double x_end, double h){
	// declare vars
	double value = 0.0;
	unsigned int nsteps = (unsigned int) ((x_end - x_start)/h);
	
	// add up values
	for (unsigned int i=1; i<nsteps; i++){
		value += h*function_handle(x_start + i*h);
	}
	value += h/2*function_handle(x_start);
	value += h/2*function_handle(x_end);
	
	return value;
}

double Integrator::simpson(double (*function_handle)(double), double x_start, double x_end, double h){
	// declare vars
	double value = 0.0;
        unsigned int ninter = (unsigned int) ((x_end - x_start)/h);
	unsigned int nsteps = ninter + 1;	

	// add up values
	for (unsigned int i=1; i<(nsteps+1)/2; i++){
                value += h/3*4*function_handle(x_start + (2*i-1)*h) + h/3*2*function_handle(x_start + 2*i*h);
        }
	value += h/3*function_handle(x_start);
	value += h/3*function_handle(x_end);

	return value;
}
