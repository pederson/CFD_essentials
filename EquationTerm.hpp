#ifndef _EQUATIONTERM_H
#define _EQUATIONTERM_H

enum EquationTermType{EQUATION_TERM_CONSTANT, EQUATION_TERM_MAIN_VARIABLE, EQUATION_TERM_VARIABLE, EQUATION_TERM_DIFFERENTIAL_OPERATOR}

class EquationTerm{
public:

	EquationTerm();
	~EquationTerm();

	void print_summary();

private:
	EquationTermType _termtype;

};

#endif