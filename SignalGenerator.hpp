#ifndef _SIGNALGENERATOR_H
#define _SIGNALGENERATOR_H

#include <math.h>

const double SIGNALGENERATOR_PI=3.14159265358979323846264338327950288;

enum class Signal : char{SIGNAL_GAUSSIAN, 
						 SIGNAL_SINUSOID, 
						 SIGNAL_CONSTANT,
						 SIGNAL_RAMP,
						 SIGNAL_TANH};

class SignalGenerator{
public:
	SignalGenerator();

	// inspectors
	double value(double t) const;
	double xloc() const {return _xloc;};
	double yloc() const {return _yloc;};
	double zloc() const {return _zloc;};


	// mutators
	void set_location(double xloc, double yloc=0.0, double zloc=0.0);
	void set_magnitude(double magnitude);
	void set_gaussian(double sigma, double t0=0.0);
	void set_sinusoid(double freq, double phase=0.0);
	void set_constant(double magnitude);
	void set_ramp(double t0, double t1);
	void set_tanh();

	

protected:

private:
	Signal _type;

	// general
	double _xloc, _yloc, _zloc;
	double _magnitude;

	// gaussian specific
	double _sigma, _t0;

	// sinusoid specific
	double _freq_Hz, _phase_deg;

	// constant specific

	// ramp specific
	double _t1;

	// tanh specific

};

#endif