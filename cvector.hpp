#ifndef _CVECTOR_H
#define _CVECTOR_H

#include <vector>
#include <iostream>

// "Composite Vector" Hey Murphy
// this is a container that holds a multiplication or division of different arrays of 
// the same length. There is also an additive constant 
class cvector{
public:

	cvector() {};
	~cvector() {};

	double operator[] (unsigned int i) const {return at(i);};
	double at(unsigned int i) const {
		double retval=1.0;
		for (auto j=0; j<_m_dvec.size(); j++) retval *= _m_dvec.at(j)[i];
		for (auto j=0; j<_m_doub.size(); j++) retval *= _m_doub.at(j);
		for (auto j=0; j<_d_dvec.size(); j++) retval *= _d_dvec.at(j)[i];
		for (auto j=0; j<_d_doub.size(); j++) retval *= _d_doub.at(j);
		return retval+_additive_const;
	}

	void multiply(const double * vec) {_m_dvec.push_back(vec);};
	//void multiply(cvector cvec);
	void multiply(const double & mval) {_m_doub.push_back(mval);};
	void divide(const double * vec) {_d_dvec.push_back(vec);};
	//void divide(cvector cvec);
	void divide(const double & dval) {_d_doub.push_back(dval);};
	void set_additive_constant(double addval) {_additive_const = addval;};

private:
	double _additive_const;

	std::vector<const double *> _m_dvec;
	//std::vector<cvector> _m_cvec;
	std::vector<double> _m_doub;
	std::vector<const double *> _d_dvec;
	std::vector<double> _d_doub;


};


#endif