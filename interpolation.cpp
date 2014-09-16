#include "interpolation.hpp"

double interp_polynomial(double * points, double * values, unsigned int degree, double eval_point){
	double interp_val=0.0, term=1.0;

	// there should be degree + 1 points and values
	for (unsigned int i=0; i<degree+1; i++){
		for (unsigned int j=0; j<degree+1; j++){
			if (j == i) continue;
			term *= (eval_point - points[j])/(points[i] - points[j]);
		}
		term *= values[i];

		interp_val += term;
		term = 1.0
	}

	return interp_val;

}

double * interp_polynomial_piecewise(double * points, double * values, unsigned int degree,
							double * eval_points, unsigned int num_eval_pts){
	// declare vars
	double * interp_vals;

	// check the input to make sure it has the correct number of inpoint points and values for 
	// the requested piecewise polynomial

	// sort the input points in ascending order and carry the index

	// sort the values by the index

	// evaluate each evaluation point for the given degree... give a pointer to 
	// the first of the interpolation points

	return interp_vals;
}