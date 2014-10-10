#ifndef _DIFFERENCE_H
#define _DIFFERENCE_H

#include <stdio.h>
#include <iostream>
#include <math.h>

class Diff{
public:

  // constructor
  Diff();

  // copy constructor

  // destructor

  // data members

  // member functions
  double forward(double fj, double fjp1, double h);
  
  double backward(double fj, double fjm1, double h);

  double central(double fjm1, double fjp1, double h);

  //dVector diff_pade(dMatrix LHS, dVector RHS);

  double second_central(double fjm1, double fj, double fjp1, double h);

protected:

private:
}
#endif 
