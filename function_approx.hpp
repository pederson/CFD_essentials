#ifndef _FUNC_APPROX_H
#define _FUNC_APPROX_H

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

enum Basis{CONSTANT=0, LINEAR=1, QUADRATIC=2, CHEBYSHEV = 3, LAGRANGE_POLY=4};

class FuncApprox{
public:

  // constructor
  FuncApprox();

  // destructor
  ~FuncApprox();

  // member data
  Basis basis_type; 
  unsigned int num_basis;
  double *a_coeff;

  // member functions
  double evaluate_basis(double *xvals, unsigned int npts);

protected:

  void basis_cheby(double (*function_handle)(double), unsigned int M);

  void basis_lagrange_poly(double (*function_handle)(double), unsigned int M);

  double * cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts);


private:
};
#endif
