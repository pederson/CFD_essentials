#ifndef _LINALGWRAPPER_H
#define _LINALGWRAPPER_H

#include <petscksp.h>
#include <iostream>
#include <vector>


void clearLinAlg(){
	PetscInt state;
	PetscFinalized(&state);
	if (!state) PetscFinalize();
}

class LinVector{
public:
	// constructor, destructor, etc...	
	LinVector(unsigned int length) {
		LinVector(); 
		set_size(length); 
		initialize(0.0); 
		VecSetFromOptions(_vec);
	}
	~LinVector() { VecDestroy(_vec);};
	LinVector(const LinVector & vec);
	operator[]{};
	operator={};

	// inspectors
	unsigned int length() const { return _length;};
	const double & data() const {
		
	} ;
	bool assembled() const { return _assembled;};

	// mutators
	void initialize(double value) { VecSet(_vec, value);};
	void insert_values(unsigned int nvals, double & vals, unsigned int & inds = nullptr){
		if (inds==nullptr && nvals==length()){
			PetscInt ix[nvals];
			for (auto i=0; i<nvals; i++) ix[i] = i;
			VecSetValues(_vec, nvals, ix, &vals, INSERT_VALUES);
		}
		else if (nvals == length()){
			VecSetValues(_vec, nvals, inds, &vals, INSERT_VALUES);
		}
		else {
			std::cout << "WARNING: indices are not specified for LinVector... values will not be set!" << std::endl;
		}
	}
	void add_values(unsigned int nvals, double & vals, unsigned int & inds = nullptr);
	void assemble() {VecAssemblyBegin(_vec); VecAssemblyEnd(_vec); _assembled = true;};



	LinVector duplicate(){
		LinVector outvec;
		VecDuplicate(_vec, &outvec._vec);
		return outvec;
	}

	LinVector subset(unsigned int start_ind, unsigned int end_ind){
		LinVector outvec;
		//im not sure how to implement this best
		return outvec;
	}

	void print_summary() {VecView(_vec,PETSC_VIEWER_STDOUT_WORLD);};

protected:

	LinVector(){
		int mpi_init;
		MPI_Initialized(&mpi_init);
		if (!mpi_init) MPI_Init(NULL, NULL);

		PetscBool init;
		PetscInitialized(&init);
		if (!init) PetscInitializeNoArguments();

		VecCreate(PETSC_COMM_WORLD,&_vec);
	}
	
	void set_size(unsigned int length){
		VecSetSizes(_vec, PETSC_DECIDE, length);
		_data.assign(length, 0.0);
		_length = length;
	}

private:

	bool _assembled;
	unsigned int _length;

	Vec _vec; // petsc vector type
	std::vector<double> _data;

};

class LinMatrix{
public:
	// constructor, destructor, etc...
	LinMatrix(unsigned int rows, unsigned int cols) {LinMatrix(); set_sizes(rows, cols); MatSetFromOptions(_mat);};
	~LinMatrix(){
		MatDestroy(&_mat);
	}
	LinMatrix(const LinMatrix & mat);

	
	// inspectors
	unsigned int nrows() const { return _nrows;};
	unsigned int ncols() const { return _ncols;};
	bool assembled() const { return _assembled;};
	void draw() const {
		MatView(_mat, PETSC_VIEWER_DRAW_WORLD);
		std::cout << "Enter anything to break the viewing: " << std::endl;
    	string input = "";
    	getline(std::cin, input);
    }

	// mutators
	void set_values(unsigned int nvals, double & vals, unsigned int & inds = NULL){

	}
	void insert_stencil();
	void add_stencil();
	void set_values_row_major(double & vals);
	void set_values_col_major(double & vals);
	void get_values();
	void assemble() {MatAssemblyBegin(_mat); MatAssemblyEnd(_mat);};
	

	set_row();
	set_column();
	get_row();
	get_column();

protected:
	LinMatrix(){
		int mpi_init;
		MPI_Initialized(&mpi_init);
		if (!mpi_init) MPI_Init(NULL, NULL);

		PetscBool init;
		PetscInitialized(&init);
		if (!init) PetscInitializeNoArguments();

		MatCreate(PETSC_COMM_WORLD,&_mat);
	}

	void set_sizes(unsigned int rows, unsigned int cols){
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
		_nrows = rows;
		_ncols = cols;
	}

	

private:

	bool _assembled;
	unsigned int _nrows, _ncols;
	Mat _mat; // petsc matrix type

};

class LinSolver{
public:
	LinSolver(){ KSPCreate(PETSC_COMM_WORLD, &_ksp);};
	~LinSolver();

	solve(LinMatrix mat, LinVector vec);

private:
	const LinMatrix & _mat;
	const LinMatrix & _pc;

	KSP _ksp; // petsc linear solver type

};

class NonlinSolver{
public:

private:

};

/*
class LinPreconditioner{
public:
	LinPreconditioner();
	~LinPreconditioner();

private:

};
*/

#endif