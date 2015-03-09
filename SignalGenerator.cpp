#include "SignalGenerator.hpp"

using namespace std;
//using namespace Signal;

SignalGenerator::SignalGenerator(){
	// set defaults
	_type = Signal::SIGNAL_CONSTANT;
	_xloc = 0.0;
	_yloc = 0.0;
	_zloc = 0.0;

	_magnitude = 1.0;

	// gaussian specific
	_sigma=0.0;
	_t0=0.0;

	// sinusoid specific
	_freq_Hz=0.0;
	_phase_deg=0.0;

	// constant specific

	// ramp specific
	_t1=0.0;
}



// mutators
void SignalGenerator::set_location(double xloc, double yloc, double zloc){
	_xloc = xloc;
	_yloc = yloc;
	_zloc = zloc;
}

void SignalGenerator::set_magnitude(double magnitude){
	_magnitude = magnitude;
}

void SignalGenerator::set_gaussian(double sigma, double t0){
	_sigma = sigma;
	_t0 = t0;
	_type = Signal::SIGNAL_GAUSSIAN;
}

void SignalGenerator::set_sinusoid(double freq, double phase){
	_freq_Hz = freq;
	_phase_deg = phase;
	_type = Signal::SIGNAL_SINUSOID;
}

void SignalGenerator::set_constant(double magnitude){
	_magnitude = magnitude;
	_type = Signal::SIGNAL_CONSTANT;
}

void SignalGenerator::set_ramp(double t0, double t1){
	_t0 = t0;
	_t1 = t1;
	_type = Signal::SIGNAL_RAMP;
}

void SignalGenerator::set_tanh(){

}

double SignalGenerator::value(double t) const{

	switch (_type){
		case Signal::SIGNAL_CONSTANT:
			return _magnitude;

		case Signal::SIGNAL_GAUSSIAN:
			return _magnitude*exp(-0.5*(t-_t0)*(t-_t0)/_sigma/_sigma);

		case Signal::SIGNAL_SINUSOID:
			return _magnitude*sin(2*SIGNALGENERATOR_PI*_freq_Hz*t + _phase_deg);

		case Signal::SIGNAL_RAMP:
			if (t>_t1) return _magnitude;
			else if (t<_t0) return 0.0;
			else return _magnitude/(_t1-_t0)*t;

		case Signal::SIGNAL_TANH:
			return _magnitude*tanh(t);

		otherwise:
			return 0.0;
	}
	
}
