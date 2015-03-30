#ifndef _SIMULATIONDATA_H
#define _SIMULATIONDATA_H

#include "Mesh.hpp"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <algorithm>

#include <H5Cpp.h>
//#include <H5File.h>
//class Static_Mesh; // forward declaration


// all data fields over a single time snapshot
class DataSnapshot{
public:
	DataSnapshot();
	~DataSnapshot();

	const double & field(std::string fieldname) const {return _datafields.at(fieldname).front();};

	friend class SimulationData;

private:
	//std::vector<std::string> _fieldnames;
	std::map<std::string, std::vector<double>> _datafields;


};

// data over the whole simulation time domain
class SimulationData{
public:

	SimulationData();
	~SimulationData();

	// inspectors
	void print_summary() const;
	bool field_present(std::string fieldname) const;
	unsigned int num_time_steps() const {return _time.size();};
	const double & get_data_at_time(double t, std::string fieldname) const {for (auto i=0; i<_time.size(); i++){ if (_time.at(i) == t){return _datasnapshots.at(i).field(fieldname);}}};
	const double & get_data_at_index(unsigned int t_index, std::string fieldname) const {return _datasnapshots.at(t_index).field(fieldname);};
	const DataSnapshot & snapshot(unsigned int t_index) const {return _datasnapshots.at(t_index);};
	const Mesh * mesh() const {return _mesh;};
	std::vector<std::string> fieldnames() const {return _fieldnames;};
	
	// mutators
	void set_time_span(double tstart, double dt, double tstop);
	void bind_mesh(const Mesh & mesh);
	void add_field(std::string fieldname);
	void add_data_at_time(double time, std::string fieldname, const double & values);
	void add_data_at_index(unsigned int t_index, std::string fieldname, const double & values);

	void write_HDF5(std::string outname) const;
	static SimulationData read_HDF5(std::string filename);

	static SimulationData combine(std::vector<const SimulationData *> datavec);

private:
	void allocate_snapshots_mesh();
	void allocate_snapshots_time();
	void allocate_snapshots_field(std::string fieldname);

	// metadata
	bool _time_set, _mesh_set;
	double _tstart, _dt, _tstop;

	// data
	const Mesh * _mesh;
	std::vector<std::string> _fieldnames;
	std::vector<DataSnapshot> _datasnapshots;

	// time info
	std::vector<double> _time;

	// used only for HDF5 reading
	Mesh _internalmesh;

	void read_HDF5_internal(std::string filename);
	
};


//H5::H5G_iterate_t get_dataset_names_H5(H5::hid_t loc_id, const char * name, const H5::H5L_info_t *linfo, void *opdata);
herr_t get_dataset_names_H5(hid_t loc_id, const char * name, const H5L_info_t *linfo, void *opdata);


#endif