#ifndef _DIFFERENCE_H
#define _DIFFERENCE_H

#include <stdio.h>
#include <iostream>
#include <math.h>

double diff_forward(double fj, double fjp1, double h);

double diff_backward(double fj, double fjm1, double h);

double diff_central(double fjm1, double fjp1, double h);

//dVector diff_pade(dMatrix LHS, dVector RHS);

double diff2_central(double fjm1, double fj, double fjp1, double h);

#endif 
