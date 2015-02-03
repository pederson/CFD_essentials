#ifndef _SIMULATIONDATA_H
#define _SIMULATIONDATA_H

#include "Mesh.hpp"

#include <vector>
#include <string>
#include <map>

/*
// data for a single time snapshot
class DataField{
public:
	DataField();
	~DataField();

	double data(i);

private:
	std::vector<double> _data;

};
*/

// all data fields over a single time snapshot
class DataSnapshot{
public:
	DataSnapshot();
	~DataSnapshot();

private:
	//std::vector<std::string> _fieldnames;
	std::map<std::string, std::vector<double>> _datafields;


};

// data over the whole simulation time domain
class SimulationData{
public:

	SimulationData();
	~SimulationData();

	DataSnapshot & at_time(double t) const {for (auto i=0; i<_time.size(); i++){ if (_time.at(i) == t){return _datasnapshots.at(i);}}};
	DataSnapshot & snapshot(unsigned int t_index) const {return _datasnapshots.at(t_index);};
	bool field_present(std::string fieldname);
	
	set_time_span(double tstart, double dt, double tstop);
	bind_mesh(const Static_Mesh & mesh);
	add_field(std::string fieldname);
	add_data(std::string fieldname, const double & values);

private:
	// metadata
	double _tstart, _dt, _tstop;
	const Static_Mesh * _mesh;
	std::vector<std::string> _fieldnames;
	//std::map<std::string, std::vector<DataSnapshot> _datasnapshots;
	std::vector<DataSnapshot> _datasnapshots;

	// time info
	std::vector<double> _time;
};


#endif