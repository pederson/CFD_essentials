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
enum PreDefinedEquation{EQUATION_ELECTROSTATIC_POISSON, 
						EQUATION_FLUID_CONTINUITY,
						EQUATION_FLUID_MOMENTUM,
						EQUATION_FLUID_ENERGY,
						EQUATION_FLUID_CONVECTIVE_DIFFUSION, 
						EQUATION_FLUID_DIFFUSION,
						EQUATION_PLASMA_CONTINUITY};

class Equation{
public:

	Equation();
	Equation(PreDefinedEquation eq);

	// inspectors
	void print_summary() const;

	// mutators
	void add_term_lhs(EquationTerm lhst);
	void add_term_rhs(EquationTerm rhst);

private:
	std::vector<EquationTerm> _lhs, _rhs;

};

#endif
