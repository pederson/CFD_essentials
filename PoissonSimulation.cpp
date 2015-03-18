#include "PoissonSimulation.hpp"

using namespace std;


PoissonSimulation::PoissonSimulation(){
	_mesh = nullptr;
	_rhs = nullptr;
	_dx = 0.0;
	_is_allocated = false;
}

PoissonSimulation::~PoissonSimulation(){

}


// mutators
//void set_boundary(BoundaryLocation loc, BoundaryCondition type, unsigned int num_layers=1);
void PoissonSimulation::bind_mesh(const RegularMesh & mesh){
	_mesh = &mesh;
	_dx = mesh.res();
}

void PoissonSimulation::bind_rhs(const cvector & rhs){
	_rhs = &rhs;
}

// other
void PoissonSimulation::view_results(){

	mesh_visualixer mvis;
	mvis.bind_mesh(*_mesh);
	mvis.set_color_ramp(CRamp::MATLAB_PARULA);
	mvis.set_colorby(&_potential.front());
	mvis.run();

}

void PoissonSimulation::run(){

	if (!_is_allocated) {
		preRunCheck();
		allocate_fields();
		_is_allocated = true;
	}

	if (_mesh->num_dims() == 1){
		cout << "1D Poisson Solver is not yet plugged in!" << endl;
	}
	else if (_mesh->num_dims() == 2){
		run_2D();
	}
	else if (_mesh->num_dims() == 3){
		cout << "3D Poisson Solver is not yet implemented!" << endl;
	}

}


void PoissonSimulation::preRunCheck(){

}

void PoissonSimulation::allocate_fields(){

	_potential.assign(_mesh->nodecount(), 0.0);

	if (_rhs == nullptr){
		_default_rhs.multiply(0.0);
		_rhs = &_default_rhs;
	}

}

void PoissonSimulation::run_2D(){

	LinVector _rhsv(_mesh->nodecount());   // approx solution, rhs 
    LinMatrix _A(_mesh->nodecount(), _mesh->nodecount());     // linear system matrix 
    LinSolver _krylov;   // linear solver context 
    int size;

    

    //cout << "about to init matrix" << endl;

	// set the matrix values for the laplace operator
	unsigned int cind, lind, rind, uind, dind;
	double oper[5];
	oper[0] = -1.0; oper[1] = -1.0; oper[2] = 4.0; oper[3] = -1.0; oper[4] = -1.0;
	int cols[5], row;
	// for each row
	for (auto j=1; j<_mesh->reg_num_nodes_x()-1; j++){ // cols
		for (auto i=1; i<_mesh->reg_num_nodes_y()-1; i++){ // rows
			cind = _mesh->reg_inds_to_glob_ind(j, i);
			lind = _mesh->reg_inds_to_glob_ind(j, i-1);
			rind = _mesh->reg_inds_to_glob_ind(j, i+1);
			uind = _mesh->reg_inds_to_glob_ind(j+1, i);
			dind = _mesh->reg_inds_to_glob_ind(j-1, i);
			cols[0] = lind;
			cols[1] = dind;
			cols[2] = cind;
			cols[3] = uind;
			cols[4] = rind;

			//if (j%10 == 0) cout << "on row: " << j << " / " << _mesh->reg_num_nodes_x() << " \r" << flush;

			row = cind;
			_A.insert_values(1, &row, 5, cols, oper);
		}
	}
	cout << endl;

	double rhv;
	for(auto i=0; i<_mesh->nodecount(); i++){
		rhv = _rhs->at(i);
		_rhsv.insert_values(1,&rhv, &i);
	}
    //cout << "assembled rhs" << endl;
	
	LinVector _soln = _krylov.solve(_A, _rhsv);
    //cout << "solved equation" << endl;
    const double * dat;
    dat = &_soln.data();
    /*
    for (auto i=0; i<100; i++){
    	cout << "rhs[" << i << "]: " << _rhsv[i] << "    " ;
    	cout << "solution[" << i << "]: " << _soln[i] << "    ";
    	cout << "solution[" << i << "]: " << dat[i] << endl;
    }
    */

    for (auto i=0; i<_mesh->nodecount(); i++) _potential[i] = dat[i];

}

