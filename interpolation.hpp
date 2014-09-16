#ifndef INTERPOLATION_H
#define INTERPOLATION_H

double interp_polynomial(double * points, double * values, unsigned int degree, double eval_point);

double * interp_polynomial_piecewise(double * points, double * values, unsigned int degree,
							double * eval_points, unsigned int num_eval_pts);

#endif