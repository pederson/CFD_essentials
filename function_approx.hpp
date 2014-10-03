#ifndef _FUNC_APPROX_H
#define _FUNC_APPROX_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

double * cheby_basis(double (*function_handle)(double), unsigned int M);

double * approx_cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts);

#endif
