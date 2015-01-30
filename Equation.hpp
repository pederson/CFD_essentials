#ifndef _EQUATION_H
#define _EQUATION_H

// The equation object should be made up of terms on the lhs and
// terms on the rhs. In turn, each term is made up of components multiplied
// together or divided. These components could be constants, variables, 
// or operators

// In addition, I should be able to define some standard equations (like
// poisson's equation) to set that equation easily

// Also, "derived" equations that calculate a quantity given other variables
// from other equations should be possible
enum PreDefinedEquation{EQUATION_POISSON, EQUATION_CONVECTIVE_DIFFUSION}
class Equation{
public:

	Equation(PreDefinedEquation eq);
	~Equation();

	void print_summary();

private:
	std::vector<EquationTerm> _lhs, _rhs;

};

#endif
