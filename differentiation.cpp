
#include "differentiation.hpp"

Diff::Diff(){
	
}

double Diff::forward(double fj, double fjp1, double h) return (fjp1 - fj)/h;

double Diff::backward(double fj, double fjm1, double h) return (fj - fjm1)/h;

double Diff::central(double fjm1, double fjp1, double h) return (fjp1 - fjm1)/(2*h);

double Diff::second_central(double fjm1, double fj, double fjp1, double h){
	return (fjp1 - 2*fj + fjm1)/(h*h);
}
