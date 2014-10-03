#define _TEST_

#define PI 3.14159265359

double * cheby_basis(double (*function_handle)(double), unsigned int M){
	// declare vars
	unsigned int *ki;
	double * si, *ak;
	double sgn;

	// initialize all s values
	si = new double[M];
	for (unsigned int i=0; i<M; i++) si[i] = i*PI/M;

	// initialize all k values
	ki = new unsigned int[M];
	for (unsigned int i=0; i<M; i++) ki[i] = i;

	// evaluate the ak coefficients by the trapezoid rule
	ak = new double[M];
	for (unsigned int i=0; i<M; i++) ak[i] = 0.0;
	for (unsigned int i=0; i<M; i++){ // loop over ki
		// determine the sign of the q(pi)=f(1)
		if (ki[i]%2 == 0) sgn = 1.0;
		else sgn = -1.0;
		ak[i] = -(function_handle(-1.0) + sgn*function_handle(1.0))/2.0;
		for (unsigned int j=; j<M; j++){ // loop over si
			ak[i] += function_handle(-cos(si[i]))*cos(ki[i]*si[i]);
		}
		ak[i] *= 2.0/double(M);
	}

	// delete stuff
	delete[] si;
	delete[] ki;

	return ak;
}

double * approx_cheby(double (*function_handle)(double), unsigned int M, double * eval_points, unsigned int N_eval_pts){
	// declare vars
	unsigned int *ki;
	double *ak, *values;

	// initialize stuff
	values = new double[N_eval_points];
	ki = new unsigned int[M];
	for (unsigned int i=0; i<M; i++) ki[i] = i;

	// get the ak coefficients
	ak = cheby_basis(function_handle, M);

	// evaluate input points
	for (unsigned int i=0; i<N_eval_points; i++){
		values[i] = ak[0]/2.0;
		for (unsigned int j=1; j<M; j++){
			values[i] += ak[j]*cos(ki[j]*acos(-eval_points[i]));
		}
	}			

	//delete stuff
	delete[] ak;
	delete[] ki;

	return values;
}

#ifdef _TEST_


#endif
