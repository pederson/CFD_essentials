#ifndef _INTEGRATION_H
#define _INTEGRATION_H

#include <math.h>
#include <stdlib.h>
#include <iostream>

class Integrator{
public:

  //constructor
  Integrator();

  // copy constructor


  // destructor
  ~Integrator();

double trapezoid(double (*function_handle)(double), double x_start, double x_end, double h);

double simpson(double (*function_handle)(double), double x_start, double x_end, double h);


}
#endif
