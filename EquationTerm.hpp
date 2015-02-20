#ifndef _EQUATIONTERM_H
#define _EQUATIONTERM_H

#include <string>
#include <vector>
#include <iostream>

enum EquationVarType{EQUATION_VAR_CONSTANT,
					  EQUATION_VAR_MAIN_VARIABLE, 
					  EQUATION_VAR_VARIABLE, 
					  EQUATION_VAR_DIFFERENTIAL_X,
					  EQUATION_VAR_DIFFERENTIAL_Y,
					  EQUATION_VAR_DIFFERENTIAL_Z,
					  EQUATION_VAR_DIFFERENTIAL_T};

enum EquationTermType{EQUATION_TERM_MULTIPLY,
					  EQUATION_TERM_DIVIDE};

class EquationTerm;

class EquationVar{
public:
	EquationVar(EquationVarType type, std::string name, double val);
	EquationVar(EquationVarType type, std::string name);
	EquationVar(EquationVarType type);


	friend class EquationTerm;

private:
	EquationVarType _vartype;
	std::string _varname;

	// optional
	double _value; // here for constants
};


class EquationTerm{
public:

	// inspectors
	void print_summary() const;

	// mutators
	void add_var(EquationVar var, EquationTermType type);

private:
	std::vector<EquationTermType> _types;
	std::vector<EquationVar> _vars;

};



#endif