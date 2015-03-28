#ifndef _POISSONSIMULATION_H
#define _POISSONSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerMesh.hpp"
#include "LinAlgWrapper.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <functional>

class PoissonSimulation{
public:

	PoissonSimulation();
	~PoissonSimulation();

	// inspectors
	const double & potential() const {return _potential.front();};


	// mutators
	void bind_mesh(const RegularMesh & mesh);
	void bind_rhs(const double * rhs);
	void bind_rhs(std::function<double(unsigned int)> rhs_fn);

	// other
	void view_results();

	void run();

private:

	void preRunCheck();
	void allocate_fields();

	void run_2D();

	const RegularMesh * _mesh;

	// user defined data
	const double * _rhs;


	// internal simulation variables
		bool _is_allocated;
		double _dx;

		// potential field
		std::vector<double> _potential;	

		std::function<double(unsigned int)> _rhs_fn;

};


#endif
