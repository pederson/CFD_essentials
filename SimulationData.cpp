#include "SimulationData.hpp"

using namespace std;

DataField::DataField();
DataField::~DataField();

double DataField::data(i);

DataSnapshot::DataSnapshot();
DataSnapshot::~DataSnapshot();

SimulationData::SimulationData();
SimulationData::~SimulationData();

bool SimulationData::field_present(std::string fieldname){
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
	_tstop = t;

	return;
}

void SimulationData::bind_mesh(const Static_Mesh & mesh){
	_mesh = &mesh; // does this copy? or just reference?
}

void SimulationData::add_field(std::string fieldname){
		if (field_present(fieldname)){
			cout << "WARNING: That field already exists!" << endl;
			return;
		}

	_fieldnames.push_back(fieldname);
	return;
}

void SimulationData::add_data(std::string fieldname, const double & values){
	if (!field_present(fieldname)){
		cout << "ERROR: that field doesn't exist" << endl;
		throw -1;
	}

	// find the correct time snapshot to insert the data
	_snapshots.at()

	for (auto i=0; i<_mesh->nodecount(); i++){

	}
}

