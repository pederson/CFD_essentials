#include "interpolation.hpp"

using namespace std;

#define _TEST_INTERPOLATE_

template <typename T1, typename T2>
void quicksort(T1 * tosort, T2 * carry, unsigned int left, unsigned int right){
	T1 pivot, tmp;
	T2 tmp2;
	unsigned int i=left, j=right;

	pivot = tosort[(right - left)/2 + left + 1];

	// partition part
	if (carry != NULL){
		while (i <= j){
			while (tosort[i] < pivot) i++;
			while (tosort[j] > pivot) j--;
			if (i <= j){
				tmp = tosort[i];
				tosort[i] = tosort[j];
				tosort[j] = tmp;
				tmp2 = carry[i];
				carry[i] = carry[j];
				carry[j] = tmp2;
				i++;
				j--;
			}
		}
	}
	else{
		while (i <= j){
			while (tosort[i] < pivot) i++;
			while (tosort[j] > pivot) j--;
			if (i <= j){
				tmp = tosort[i];
				tosort[i] = tosort[j];
				tosort[j] = tmp;
			
				i++;
				j--;
			}
		}
	}

	// recursion part
	if (left < j) quicksort(tosort, carry, left, j);
	if (i < right) quicksort(tosort, carry, i, right);

	return;
}

double interp_nearest(double * points, double * values, unsigned int numpts, double eval_point){
	// This assumes that the points are sorted in increasing order

	// declare vars
	double min_diff, ldiff, rdiff;
	unsigned int i;
	bool cont = true;

	// find a starting point
	i = numpts * (eval_point - points[0])/(points[numpts-1] - points[0]);
	min_diff = (points[i] - eval_point)*(points[i] - eval_point);

cout << "Starting point: " << points[i] << " mindiff: " << min_diff << endl;

	// iterate to find the closest point
	while (cont){
		// check the current point for the bounds
		if (i == 0){
			if ((points[1] - eval_point)*(points[1] - eval_point) < min_diff){
				i = 1;
				min_diff = (points[1] - eval_point)*(points[1] - eval_point);
			}
			else break;
		}
		if (i == numpts-1){
			if ((points[numpts-2] - eval_point)*(points[numpts-2] - eval_point) < min_diff){
				i = numpts-2;
				min_diff = (points[numpts-2] - eval_point)*(points[numpts-2] - eval_point);
			}
			else break;
		}

		// check the current point left and right and step accordingly
		ldiff = (points[i-1] - eval_point)*(points[i-1] - eval_point);
		rdiff = (points[i+1] - eval_point)*(points[i+1] - eval_point);
		if (ldiff < min_diff){
			min_diff = ldiff;
			i--;
			continue;
		}
		if (rdiff < min_diff){
			min_diff = rdiff;
			i++;
			continue;
		}

		break;

	}

	return values[i];
}

double * interp_nearest_piecewise(double * points, double * values, double * eval_points,
							unsigned int num_eval_pts){
	// declare vars
	double * interp_vals;

	// check the input points

	// sort the input points in ascending order and carry the index

	// sort the values by the index

	// evaluate the nearest neighbor for each eval point


	return interp_vals;
}

double interp_polynomial(double * points, double * values, unsigned int degree, double eval_point){
	double interp_val=0.0, term=1.0;

	// there should be degree + 1 points and values
	for (unsigned int i=0; i<degree+1; i++){
		for (unsigned int j=0; j<degree+1; j++){
			if (j == i) continue;
			term *= (eval_point - points[j])/(points[i] - points[j]);
		}
		term *= values[i];

		interp_val += term;
		term = 1.0;
	}

	return interp_val;

}

double * interp_polynomial_piecewise(double * points, double * values, unsigned int degree,
							double * eval_points, unsigned int num_eval_pts){
	// declare vars
	double * interp_vals;

	// check the input to make sure it has the correct number of inpoint points and values for 
	// the requested piecewise polynomial

	// sort the input points in ascending order and carry the index

	// sort the values by the index

	// evaluate each evaluation point for the given degree... give a pointer to 
	// the first of the interpolation points

	return interp_vals;
}

#ifdef _TEST_INTERPOLATE_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char * argv[]){
	// vars
	double * randpts, * randvals, * tosort;
	double evalpt, nnval, linval, quadval, cubval;
	int npts;

	// user set variables
	npts = 50;

	// generate randpts and randvals;
	srand(time(NULL));
	randpts = new double[npts];
	randvals = new double[npts];
	for (int i=0; i<npts; i++){
		randpts[i] = rand()%npts + 1;
		randvals[i] = rand()%20 - 10;
	}

	// test quicksort
	cout << "testing quicksort" << endl;
	tosort = new double[npts];
	for (int i=0; i<npts; i++) tosort[i] = randpts[i];
	quicksort(tosort, (int *) NULL, 0, npts-1);
	for (int i=0; i<npts; i++){
		cout << "unsorted: " << randpts[i] << " sorted: " << tosort[i] << endl;
	}
	quicksort(randpts, randvals, 0, npts-1);

	// test nearest neighbor interpolation
	cout << "testing nearest neighbor interpolation" << endl;
	evalpt = rand()%npts + 1 + 0.4;
	nnval = interp_nearest(randpts, randvals, npts, evalpt);
	cout << "evalpt: " << evalpt << " nearest val: " << nnval << endl;
	for (int i=0; i<npts; i++){
		cout << "pts: " << randpts[i] << " vals: " << randvals[i] << endl;
	}

	return 0;
}

#endif