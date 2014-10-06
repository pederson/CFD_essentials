
#include "differentiation.hpp"

double diff_forward(double fj, double fjp1, double h) return (fjp1 - fj)/h;

double diff_backward(double fj, double fjm1, double h) return (fj - fjm1)/h;

double diff_central(double fjm1, double fjp1, double h) return (fjp1 - fjm1)/(2*h);

double diff2_central(double fjm1, double fj, double fjp1, double h){
	return (fjp1 - 2*fj + fjm1)/(h*h);
}
