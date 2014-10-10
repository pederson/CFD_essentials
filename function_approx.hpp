#ifndef _FUNC_APPROX_H
#define _FUNC_APPROX_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

class FuncApprox{
public:

  // constructor
  FuncApprox();

  // copy constructor

  // destructor
  ~FuncApprox();

  double * basis_cheby(double (*function_handle)(double), unsigned int M);

  double * basis_lagrange_poly(double (*function_handle)(double), unsigned int M);

  double * cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts);

protected:

private:
}
#endif
