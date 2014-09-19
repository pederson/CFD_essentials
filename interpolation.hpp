#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <stdlib.h>
#include <iostream>
#include <math.h>

//template <typename T1, typename T2> quicksort(T1 tosort, T2 carry, unsigned int left, unsigned int right);

double interp_nearest(double * points, double * values, unsigned int numpts, double eval_point);

double * interp_nearest_piecewise(double * points, double * values, unsigned int npts, 
							double * eval_points, unsigned int num_eval_pts);

double interp_polynomial(double * points, double * values, unsigned int degree, double eval_point);

double * interp_polynomial_piecewise(double * points, double * values, unsigned int npts, unsigned int degree,
							double * eval_points, unsigned int num_eval_pts);

double * chebyshev(unsigned int N, double min_val=-1.0, double max_val=1.0);
#endif