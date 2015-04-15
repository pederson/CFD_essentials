// The equation object should be made up of terms on the lhs and
// terms on the rhs. In turn, each term is made up of components multiplied
// together or divided. These components could be constants, variables, 
// or operators

// In addition, I should be able to define some standard equations (like
// poisson's equation) to set that equation easily

//#define _TEST_

// Also, "derived" equations that calculate a quantity given other variables
// from other equations should be possible
enum PreDefinedEquation{EQUATION_ELECTROSTATIC_POISSON, 
						EQUATION_FLUID_CONTINUITY,
						EQUATION_FLUID_MOMENTUM,
						EQUATION_FLUID_ENERGY,
						EQUATION_FLUID_CONVECTIVE_DIFFUSION, 
						EQUATION_FLUID_DIFFUSION,
						EQUATION_PLASMA_CONTINUITY};


using namespace std;

Equation::Equation(){

}

Equation::Equation(PreDefinedEquation eq){

	switch (eq){
		case EQUATION_ELECTROSTATIC_POISSON:
			EquationTerm term1;

			term1.add_var(EquationVar(EQUATION_VAR_CONSTANT, "eps0", 1.1e-19), EQUATION_TERM_MULTIPLY);
			term1.add_var(EquationVar(EQUATION_VAR_DIFFERENTIAL_T), EQUATION_TERM_MULTIPLY);
			term1.add_var(EquationVar(EQUATION_VAR_MAIN_VARIABLE, "n_e"), EQUATION_TERM_MULTIPLY);
			term1.add_var(EquationVar(EQUATION_VAR_VARIABLE, "red"), EQUATION_TERM_DIVIDE);

			break;

		case EQUATION_FLUID_CONTINUITY:

			break;

		case EQUATION_FLUID_MOMENTUM:

			break;

		case EQUATION_FLUID_ENERGY:

			break;

		case EQUATION_FLUID_DIFFUSION:

			break;

		case EQUATION_PLASMA_CONTINUITY:

			break;

		otherwise:
			cout << "ERROR: Unrecognized predefined equation" << endl;

	}

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


#ifdef _TEST_

int main(int argc, char * argv[]){



	return 0;
}