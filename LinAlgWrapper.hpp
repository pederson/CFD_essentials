#ifndef _LINALGWRAPPER_H
#define _LINALGWRAPPER_H

#include <petscksp.h>
#include <iostream>
#include <vector>
#include <string>

//using PreCondType = PCType;
//using PRECOND_NONE = PCNONE;
//using PRECOND_JACOBI = PCJACOBI;

//enum PreCondType{PRECOND_NONE, PRECOND_JACOBI};

void clearLinAlg(){
	PetscBool state;
	PetscFinalized(&state);
	if (!state) PetscFinalize();
}

class LinVector{
public:
	// constructor, destructor, etc...	
	LinVector(unsigned int length) {
		//PetscErrorCode ierr;
		LinVector(); 
		//VecSetSizes(this->_vec, PETSC_DECIDE, length);
		set_size(length); 
		//VecSetType(_vec, VECSEQ);
		VecSetFromOptions(_vec);
		//std::cout << "error code: " << ierr << std::endl;
		
		initialize(0.0); 
		
	}
	~LinVector() { VecDestroy(&_vec);};
	//LinVector(const LinVector & vec);
	//const double & operator[]{};
	//LinVector operator={};

	// inspectors
	unsigned int length() const { return _length;};
	const double & data() const { return _data.front();};

	// this doesn't really belong in public
	bool assembled() const { return _assembled;};

	// mutators
	void initialize(double value) { VecSet(_vec, value); disassemble();};
	void insert_values(unsigned int nvals, const double * vals, const int * inds = nullptr){
		if (inds==nullptr && nvals==length()){
			PetscInt ix[nvals];
			for (auto i=0; i<nvals; i++) ix[i] = i;
			VecSetValues(_vec, nvals, ix, vals, INSERT_VALUES);
		}
		else if (nvals == length()){
			VecSetValues(_vec, nvals, inds, vals, INSERT_VALUES);
		}
		else {
			std::cout << "WARNING: indices are not specified for LinVector... values will not be set!" << std::endl;
		}
		disassemble();
	}
	void add_values(unsigned int nvals, const double * vals, const int * inds = nullptr);
	


	LinVector duplicate(){
		LinVector outvec;
		VecDuplicate(_vec, &outvec._vec);
		outvec._data = _data;
		return outvec;
	}

	//LinVector subset(unsigned int start_ind, unsigned int end_ind){
	//	LinVector outvec;
		//im not sure how to implement this best
	//	return outvec;
	//}

	void print_summary() {VecView(_vec,PETSC_VIEWER_STDOUT_WORLD);};

	friend class LinSolver;

protected:

	LinVector(){
		int mpi_init;
		MPI_Initialized(&mpi_init);
		if (!mpi_init) {
			//std::cout << "init MPI" << std::endl;
			MPI_Init(NULL, NULL);
		}

		PetscBool init;
		PetscInitialized(&init);
		if (!init) {
			//std::cout << "init PETSC" << std::endl;
			PetscInitializeNoArguments();
		}

		//PetscErrorCode ierr;
		VecCreate(PETSC_COMM_WORLD, &_vec); 
	}
	
	void set_size(unsigned int length){
		_length = length;
		VecSetSizes(_vec, PETSC_DECIDE, length);
		_data.assign(length, 0.0);
		
	}


private:

	bool _assembled;
	unsigned int _length;

	Vec _vec; // petsc vector type
	std::vector<double> _data;

	void assemble() {VecAssemblyBegin(_vec); VecAssemblyEnd(_vec); _assembled = true;};
	void disassemble() {_assembled = false;};
	void refresh_data() {
		double soln[_length];
	    PetscInt inds[_length];
	    for (auto i=0; i<_length; i++) inds[i] = i;
	    VecGetValues(_vec, _length, inds, soln);
		for (auto i=0; i<_length; i++) _data[i] = soln[i];
	}

};

class LinMatrix{
public:
	// constructor, destructor, etc...
	LinMatrix(unsigned int rows, unsigned int cols) {
		LinMatrix(); 
		set_sizes(rows, cols); 
		//MatSetFromOptions(_mat);
		MatSetType(_mat, MATMAIJ);
	}
	~LinMatrix(){ MatDestroy(&_mat);};
	LinMatrix(const LinMatrix & mat);

	
	// inspectors
	unsigned int nrows() const { return _nrows;};
	unsigned int ncols() const { return _ncols;};
	void draw() {
		if (!_assembled) assemble();
		MatView(_mat, PETSC_VIEWER_DRAW_WORLD);
		std::cout << "Enter anything to break the viewing: " << std::endl;
    	std::string input = "";
    	getline(std::cin, input);
    	disassemble();
	}

	// this doesn't really belong in public
	bool assembled() const { return _assembled;};


	// mutators
	void insert_values(unsigned int numrows, const int * rows, unsigned int numcols, const int * cols, const double * vals){
		MatSetValues(_mat, numrows, rows, numcols, cols, vals, INSERT_VALUES);
	}
	//void add_values(unsigned int nvals, double & vals, unsigned int & inds = NULL);
	//void insert_stencil(unsigned int nvals, double & vals, unsigned int & inds){

	//}
	//void add_stencil();

	//void set_values_row_major(double & vals);
	//void set_values_col_major(double & vals);
	//void get_values();

	//get_row();
	//get_column();

	friend class LinSolver;

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

	void disassemble() {_assembled = false;};
	void assemble() {MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);};

};



class LinSolver{
public:
	LinSolver(){ KSPCreate(PETSC_COMM_WORLD, &_ksp); _pctype = PCJACOBI;};
	~LinSolver(){KSPDestroy(&_ksp);};

	void setPreConditioner(PCType pctype) {_pctype = pctype;};
	//void setTolerances();

	LinVector solve(LinMatrix & mat, LinVector & rhs){
		if (!mat.assembled()) mat.assemble();
		if (!rhs.assembled()) rhs.assemble();
		LinVector outvec(rhs.length());
		KSPSetOperators(_ksp, mat._mat, mat._mat);
		KSPGetPC(_ksp, &_pc);
		PCSetType(_pc, _pctype);
		KSPSetTolerances(_ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
		KSPSolve(_ksp, rhs._vec, outvec._vec);
		outvec.refresh_data();
		return outvec;
	}

private:
	PC _pc;

	PCType _pctype;
	KSP _ksp; // petsc linear solver type

};



class NonlinSolver{
public:

private:

};


#endif