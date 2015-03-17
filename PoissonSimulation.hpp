#ifndef _POISSONSIMULATION_H
#define _POISSONSIMULATION_H

#include "RegularMesh.hpp"
#include "VisualixerMesh.hpp"
#include "cvector.hpp"
#include "LinAlgWrapper.hpp"

#include <vector>
#include <string>
#include <iostream>

class PoissonSimulation{
public:

	PoissonSimulation();
	~PoissonSimulation();

	// inspectors
	const double & potential() const {return _potential.front();};
	const cvector & rhs() const {return *_rhs;};


	// mutators
	void bind_mesh(const RegularMesh & mesh);
	void bind_rhs(const cvector & rhs);

	// other
	void view_results();

	void run();

private:

	void preRunCheck();
	void allocate_fields();

	void run_2D();

	const RegularMesh * _mesh;

	// user defined data
	const cvector * _rhs;


	// internal simulation variables
		bool _is_allocated;
		double _dx;

		// potential field
		std::vector<double> _potential;	

		// defaults
		cvector _default_rhs;

};


#endif
