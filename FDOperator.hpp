#ifndef _FDOPERATOR_H
#define _FDOPERATOR_H

#include "RegularMesh.hpp"

#include <iostream>
#include <vector>
#include <string>

enum FDOperatorType{FD_DERIVATIVE_X, FD_DERIVATIVE_Y, FD_DERIVATIVE_Z, FD_DIVERGENCE, FD_LAPLACIAN};
enum BoundaryMethod{BOUND_EXTEND, BOUND_ZERO, BOUND_PARTIAL_DERIV}

class FDOperator{
public:

private:
	const RegularMesh * m_mesh;
	FDOperatorType m_op;
	BoundaryMethod m_bmeth;

	std::vector<double> m_result;
};

#endif