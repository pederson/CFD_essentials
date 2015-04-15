#include "EquationTerm.hpp"

//#define _TEST_

using namespace std;

// inspectors
void EquationTerm::print_summary() const{
	// print all the "multiply" variables first
	unsigned int mtc=0, dtc=0, ct=0;
	for (auto i=0; i<_types.size(); i++){
		if (_types[i] == EQUATION_TERM_MULTIPLY) mtc++;
		if (_types[i] == EQUATION_TERM_DIVIDE) dtc++;
	}

	if (mtc == 0) cout << "1" ;
	else{
		cout << "(" ;
		for (auto i=0; i<_types.size(); i++){
			if (_types[i] == EQUATION_TERM_MULTIPLY){
				cout << _vars[i]._varname ;
				break;
			}
		}
		if (mtc!=1){
			ct=0;
			for (auto i=0; i<_types.size(); i++){
				if (_types[i] == EQUATION_TERM_MULTIPLY){
					if (ct>0) cout << "*" << _vars[i]._varname ;
					ct++;
				}
			}
		}
		cout << ")";
	}
	if (dtc!=0){
		cout << "/(" ;
		// deal with the first variable name
		for (auto i=0; i<_types.size(); i++){
			if (_types[i] == EQUATION_TERM_DIVIDE){
				cout << _vars[i]._varname ;
				break;
			}
		}
		// deal with the other variable names
		if (dtc!=1){
			ct=0;
			for (auto i=0; i<_types.size(); i++){
				if (_types[i] == EQUATION_TERM_DIVIDE){
					if (ct>0) cout << "*" << _vars[i]._varname ;
					ct++;
				}
			}
		}
		cout << ")" ;
	}
}

// mutators
void EquationTerm::add_var(EquationVar var, EquationTermType type){
	_vars.push_back(var);
	_types.push_back(type);
	return;
}

EquationVar::EquationVar(EquationVarType type, std::string name, double val){
	_vartype = type;
	_varname = name;
	_value = val;
}

EquationVar::EquationVar(EquationVarType type, std::string name){
	_vartype = type;
	_varname = name;
}

EquationVar::EquationVar(EquationVarType type){
	_vartype = type;

	switch (type){
		case (EQUATION_VAR_DIFFERENTIAL_X):
			_varname = "d/dx";
			break;

		case (EQUATION_VAR_DIFFERENTIAL_Y):
			_varname = "d/dy";
			break;
		
		case (EQUATION_VAR_DIFFERENTIAL_Z):
			_varname = "d/dz";
			break;
			
		case (EQUATION_VAR_DIFFERENTIAL_T):
			_varname = "d/dt";
			break;
		otherwise:
			cout << "ERROR: Unrecognized variable type" << endl;
	}
}


#ifdef _TEST_

// compile with g++ -std=c++11 EquationTerm.cpp -o equationterm_test

int main(int argc, char * argv[]){
	EquationTerm term;
	term.add_var(EquationVar(EQUATION_VAR_CONSTANT, "eps0", 1.1e-19), EQUATION_TERM_MULTIPLY);
	term.add_var(EquationVar(EQUATION_VAR_DIFFERENTIAL_T), EQUATION_TERM_MULTIPLY);
	term.add_var(EquationVar(EQUATION_VAR_MAIN_VARIABLE, "n_e"), EQUATION_TERM_MULTIPLY);
	term.add_var(EquationVar(EQUATION_VAR_VARIABLE, "red"), EQUATION_TERM_DIVIDE);

	term.print_summary();
	cout << endl;

	return 0;
}

#endif