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

Equation::Equation(){

}
Equation::Equation(PreDefinedEquation eq);
~Equation(){

}

// inspectors
void Equation::print_summary() const{
	cout << "**** Equation Summary ****" << endl;
	cout << "\t\t" << endl;

	// print lhs first
	if (_lhs.size() == 0) cout << " 0" ;
	for (auto i=0; i<_lhs.size(); i++) {
		_lhs[i].print_summary();
		cout << " + " ;
	}
	cout << " = " ;
	for (auto i=0; i<_rhs.size(); i++) {
		_rhs[i].print_summary();
		cout << " + " ;
	}
	cout << endl;
	cout << "**************************" << endl;

}

// mutators
void Equation::add_term_lhs(EquationTerm lhst){
	_lhs.push_back(lhst);
}

void Equation::add_term_rhs(EquationTerm rhst){
	_rhs.push_back(rhst);
}

private:
	std::vector<EquationTerm> _lhs, _rhs;

};

#endif
