#include "SimulationData.hpp"

using namespace std;

DataSnapshot::DataSnapshot(){

}

DataSnapshot::~DataSnapshot(){

}

SimulationData::SimulationData(){
	_mesh_set = false;
	_time_set = false;
	_tstart = 0.0;
	_dt = 0.0;
	_tstop = 0.0;
	_datasnapshots.resize(1); // allocate one snapshot
	_time.resize(1);
}

SimulationData::~SimulationData(){

}

void SimulationData::print_summary() const{
	cout << " " << endl;
	cout << "******** Simulation Data Summary ******** " << endl;
	if (_datasnapshots.size() == 0){
	cout << "  Data is empty!" << endl;
	return;
	}

	if (_mesh_set) cout << "\tMesh is bound (" << _mesh->nodecount() << " nodes; " << _mesh->elementcount() << " elements)" << endl;
	cout << "\tNumber of time steps: " << _time.size() << endl;
	cout << "\tTime Start: " << _tstart << "\tTime Step: " << _dt << "\tTime Stop: " << _tstop << endl;

	cout << "\tField Names:" << endl;
	for (auto i=0; i<_fieldnames.size(); i++) cout << "\t\t" << _fieldnames.at(i) << endl;
	cout << "************************************* " << endl;
	cout << " " << endl;
}

bool SimulationData::field_present(std::string fieldname) const {
	for (auto i=0; i<_fieldnames.size(); i++){
		if (_fieldnames.at(i).compare(fieldname) == 0){
			return true;
		}
	}
	return false;
}


void SimulationData::set_time_span(double tstart, double dt, double tstop){
	_tstart = tstart;
	_dt = dt;
	_tstop = tstop;

	_time_set = true;

	// construct the time vector
	double tcur = _tstart;
	while (tcur <= _tstop){
		_time.push_back(tcur);
		tcur += dt;
	}

	// expand the data snapshot vector accordingly
	_datasnapshots.resize(_time.size());

	// allocate resources in DataSnapshots here
	allocate_snapshots_time();
	return;
}

void SimulationData::bind_mesh(const Static_Mesh & mesh){
	_mesh = &mesh; 
	_mesh_set = true;
	allocate_snapshots_mesh();
}

void SimulationData::add_field(std::string fieldname){
		if (field_present(fieldname)){
			cout << "WARNING: That field already exists!" << endl;
			return;
		}

	_fieldnames.push_back(fieldname);

	// allocate resources in DataSnapshots here
	allocate_snapshots_field(fieldname);
	return;
}


void SimulationData::add_data_at_time(double time_s, std::string fieldname, const double & values){
	bool added = false;
	for (auto i=0; i<_time.size(); i++){ 
		if (_time.at(i) == time_s){
			add_data_at_index(i, fieldname, values);
			break;
			added = true;
		}
	}
	if (!added){
		cout << "ERROR: that time slot doesn't exist" << endl;
		throw -1;
	}
}

void SimulationData::add_data_at_index(unsigned int t_index, std::string fieldname, const double & values){
	if (!field_present(fieldname)){
		cout << "ERROR: that field doesn't exist" << endl;
		throw -1;
	}

	const double * vals = &values;
	for (auto i=0; i<_mesh->nodecount(); i++){
		_datasnapshots.at(t_index)._datafields[fieldname].at(i) = vals[i];
	}
}


// run this when a new field has been added
void SimulationData::allocate_snapshots_field(string fieldname){

	if (_mesh_set){
		// set the first time snapshot by default
		_datasnapshots.at(0)._datafields[fieldname].assign(_mesh->nodecount(), 0.0);


		// set others if more time snapshots are allocated
		if (_time_set){
			for (auto j=1; j<_datasnapshots.size(); j++){
				_datasnapshots.at(j)._datafields[fieldname].assign(_mesh->nodecount(), 0.0);
			}
		}
	}

}

// run this when a mesh has been added 
void SimulationData::allocate_snapshots_mesh(){

	string fieldname;
	for (auto i=0; i<_fieldnames.size(); i++){
		fieldname = _fieldnames.at(i);

		// set the first time snapshot by default

		_datasnapshots.at(0)._datafields[fieldname].assign(_mesh->nodecount(), 0.0);


		// set others if more time snapshots are allocated
		if (_time_set){
			for (auto j=1; j<_datasnapshots.size(); j++){
				_datasnapshots.at(j)._datafields[fieldname].assign(_mesh->nodecount(), 0.0);
			}
		}
	}

}

// run this when time span has been set 
void SimulationData::allocate_snapshots_time(){

	string fieldname;
	if (_mesh_set){
		for (auto i=0; i<_fieldnames.size(); i++){
			// set the first time snapshot by default
			fieldname = _fieldnames.at(i);


			for (auto j=0; j<_datasnapshots.size(); j++){
				_datasnapshots.at(j)._datafields[fieldname].assign(_mesh->nodecount(), 0.0);
			}
		}
	}

}

