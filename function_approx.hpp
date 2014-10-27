#ifndef _FUNC_APPROX_H
#define _FUNC_APPROX_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

enum Basis{CHEBYSHEV = 0, LAGRANGE_POLY=1, CONSTANT=2, LINEAR=3, QUADRATIC=4};

class FuncApprox{
public:

  // constructor
  FuncApprox();

  // copy constructor

  // destructor
  ~FuncApprox();

  // member data
  unsigned int num_basis;
  double *a_coeff;

  // member functions
  double * basis_cheby(double (*function_handle)(double), unsigned int M);

  double * basis_lagrange_poly(double (*function_handle)(double), unsigned int M);

  double * cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts);

protected:

private:
};
#endif
