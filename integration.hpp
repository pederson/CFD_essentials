#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include <math.h>
#include <stdlib.h>
#include <iostream>

double integrate_trapezoid(double (*function_handle)(double), double x_start, double x_end, double h);

double integrate_simpson(double (*function_handle)(double), double x_start, double x_end, double h);



#endif
