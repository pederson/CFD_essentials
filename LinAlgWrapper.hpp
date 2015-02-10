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
	LinVector(unsigned int length) {LinVector(); set_sizes(length); VecSetFromOptions(_vec);};

	void set_size(unsigned int length){VecSetSizes(_vec, PETSC_DECIDE, length);};
	unsigned int get_size() { PetscInt sz; VecGetSize(_vec, &sz); return sz;};
	
	void set(double value){VecSet(_vec, value);};
	void set_values(unsigned int nvals, double & vals, unsigned int & inds = NULL){
		if (inds==NULL && nvals==get_size()){
			PetscInt ix[nvals];
			for (auto i=0; i<nvals; i++) ix[i] = i;
			VecSetValues(_vec, nvals, ix, &vals, INSERT_VALUES);
		}
	}
	void get_values();
	void assemble() {VecAssemblyBegin(_vec); VecAssemblyEnd(_vec);};

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
	~LinVector(){
		VecDestroy(_vec);
	};

	void set_type();
	void get_type();

private:
	//bool size_set;
	//bool values_set;
	//bool type_set;
	//bool assembly_set;

	Vec _vec; // petsc vector type

};

class LinMatrix{
public:
	
	LinMatrix(unsigned int rows, unsigned int cols) {LinMatrix(); set_sizes(rows, cols); MatSetFromOptions(_mat);};

	void set_size(unsigned int rows, unsigned int cols){MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);};
	unsigned int get_size() { PetscInt nrows, ncols; MatGetSize(_vec, &nrows, &ncols); return nrows*ncols;};
	
	void set_values(unsigned int nvals, double & vals, unsigned int & inds = NULL){
		if (inds==NULL && nvals==get_size()){
			PetscInt ix[nvals];
			for (auto i=0; i<nvals; i++) ix[i] = i;
			VecSetValues(_vec, nvals, ix, &vals, INSERT_VALUES);
		}
	}
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
	~LinMatrix(){
		MatDestroy(&_mat);
	}

private:

	Mat _mat; // petsc matrix type

};

class LinSolver{
public:
	LinSolver();
	~LinSolver();

private:

	KSP _ksp; // petsc linear solver type

};

class NonlinSolver{
public:

private:

};

class LinPreconditioner{
public:
	LinPreconditioner();
	~LinPreconditioner();

private:

};

#endif