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

	_time.clear();
	// construct the time vector
	double tcur = _tstart;
	while (tcur <= _tstop){
		//cout << "Tcur: " << tcur << "\t tstop: " << tstop << endl;
		_time.push_back(tcur);
		
		tcur += dt;
	}

	// expand the data snapshot vector accordingly
	_datasnapshots.resize(_time.size());

	// allocate resources in DataSnapshots here
	allocate_snapshots_time();
	return;
}

void SimulationData::bind_mesh(const Mesh & mesh){
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

void SimulationData::write_HDF5(std::string outname) const{
	//#ifdef _HDF5_H

	cout << "Writing SimulationData into HDF5 file: " << outname << " ..... this could take a while. grab some coffee" << endl;
	// create an empty HDF5 file
	H5::H5File outfile(outname, H5F_ACC_TRUNC);

	// ******* create a group for the mesh *******
	H5::Group meshgroup(outfile.createGroup("/Mesh"));

		H5::Group nodes_group(meshgroup.createGroup("/Mesh/Nodes"));
		H5::Group elements_group(meshgroup.createGroup("/Mesh/Elements"));
		H5::Group nodedata_group(meshgroup.createGroup("/Mesh/NodeData"));
		H5::Group elementdata_group(meshgroup.createGroup("/Mesh/ElementData"));

		hsize_t nnodes = _mesh->nodecount();
		hsize_t nelem = _mesh->elementcount();
		H5::DataSpace nodespace(1, &nnodes);
		H5::DataSpace elemspace(1, &nelem);

		// insert nodes
		double nodebuf[nnodes];
		for (auto i=0; i<nnodes; i++) nodebuf[i] = _mesh->node(i).x();
		H5::DataSet nodex_set = nodes_group.createDataSet("NodesX", H5::PredType::NATIVE_DOUBLE, nodespace);
		nodex_set.write(nodebuf, H5::PredType::NATIVE_DOUBLE);
		nodex_set.close();

		if (_mesh->num_dims() > 1){
			for (auto i=0; i<nnodes; i++) nodebuf[i] = _mesh->node(i).y();
			H5::DataSet nodey_set = nodes_group.createDataSet("NodesY", H5::PredType::NATIVE_DOUBLE, nodespace);
			nodey_set.write(nodebuf, H5::PredType::NATIVE_DOUBLE);
			nodey_set.close();
		}
		
		if (_mesh->num_dims() > 2){
			for (auto i=0; i<nnodes; i++) nodebuf[i] = _mesh->node(i).z();
			H5::DataSet nodez_set = nodes_group.createDataSet("NodesZ", H5::PredType::NATIVE_DOUBLE, nodespace);
			nodez_set.write(nodebuf, H5::PredType::NATIVE_DOUBLE);
			nodez_set.close();
		}

		// insert elements

		// insert nodedata

		// insert elementdata

		// mesh cleanup
		nodespace.close();
		elemspace.close();
		elementdata_group.close();
		nodedata_group.close();
		elements_group.close();
		nodes_group.close();
		meshgroup.close();
	//cout << "WROTE MESH" << endl;

	// ******* create a group for the Fields *******
	H5::Group fields_group(outfile.createGroup("/Fields"));
	//cout << "CREATED FIELDS" << endl;

		hsize_t ntime = num_time_steps();
		hsize_t fielddim[2];
		fielddim[0] = nnodes;
		fielddim[1] = ntime;
		hsize_t fielddims[2];
		fielddims[0] = nnodes;
		fielddims[1] = 1;

		H5::DataSpace fieldspace(2, fielddim);
		H5::DataSet field_set;
		H5::DataSpace memspace = H5::DataSpace(2, fielddims);

		//cout << "CREATED FIELD DATASETS" << endl;
		
		//double ** fieldbuf = new double*[nnodes];
		//for (auto i=0; i<nnodes; i++) fieldbuf[i] = new double[ntime];
		const double * fieldbuf;
		hsize_t offset[2], count[2], stride[2], block[2];
		offset[0] = 0;

		count[0]  = nnodes;
		count[1]  = 1;

		stride[0] = 1;
		stride[1] = 1;

		block[0] = 1;
		block[1] = 1;

		string cfieldname;
		for (auto i=0; i<_fieldnames.size(); i++){
			cfieldname = _fieldnames.at(i);
			//cout << "writing field: " << cfieldname << endl;

			field_set = fields_group.createDataSet(_fieldnames.at(i), H5::PredType::NATIVE_DOUBLE, fieldspace);
			//cout << "created dataset" << endl;
			for (auto j=0; j<ntime; j++){
				cout << "writing time step " << j << "/" << ntime << "\r" << flush;
				offset[1] = j;
				fieldspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);


				fieldbuf = &_datasnapshots.at(j)._datafields.at(cfieldname).front();
				field_set.write(fieldbuf, H5::PredType::NATIVE_DOUBLE, memspace, fieldspace);
				/*
				for (auto k=0; k<nnodes; k++){
					fieldbuf[k][j] = _datasnapshots.at(j)._datafields.at(cfieldname).at(k);
					//if (j==1) cout << "j: " << j << " k: " << k << " fieldbuf: " << fieldbuf[k][j] << endl;// "\r" << flush;
				}
				*/

			}
			//cout << "wrote to buffer: " << fieldbuf[0][0] << endl;
			//field_set.write(fieldbuf, H5::PredType::NATIVE_DOUBLE);
			//cout << "wrote to file" << endl;
			field_set.close();
		}
		//for (auto i=0; i<nnodes; i++) delete[] fieldbuf[i];
		//delete[] fieldbuf;

		// fields cleanup
		fieldspace.close();
		fields_group.close();

	//cout << "WROTE FIELDS" << endl;

	// ******* create a group for the Time values *******
	H5::Group time_group(outfile.createGroup("/Time"));

		H5::DataSpace timespace(1, &ntime);
		H5::DataSet time_set = time_group.createDataSet("Time", H5::PredType::NATIVE_DOUBLE, timespace);

		double timebuf[ntime];
		for (auto i=0; i<ntime; i++) timebuf[i] = _time.at(i);

		time_set.write(timebuf, H5::PredType::NATIVE_DOUBLE);
		time_set.close();

		// time cleanup
		timespace.close();
		time_group.close();

	//cout << "WROTE TIME " << endl;

	// close the file
	outfile.close();
	//#else 
	//	cout << "HDF5 library is not included...doing nothing" << endl;
	//#endif

}

SimulationData SimulationData::read_HDF5(std::string filename){
	
	SimulationData simdat;
	simdat.read_HDF5_internal(filename);
	return simdat;

}

std::vector<std::string> HDF5FieldNames;
void SimulationData::read_HDF5_internal(std::string filename){
	cout << "Reading HDF5 file: " << filename << endl;

	// open the file for reading
	H5::H5File file(filename, H5F_ACC_RDONLY);

	// Mesh reading
	H5::Group mesh_group = file.openGroup("Mesh");

		_mesh = &_internalmesh;
		_mesh_set = true;

		H5::Group nodes_group = mesh_group.openGroup("Nodes");
		H5::Group elements_group = mesh_group.openGroup("Elements");
		H5::Group nodedata_group = mesh_group.openGroup("NodeData");
		H5::Group elementdata_group = mesh_group.openGroup("ElementData");

		H5::DataSet nodes_set = nodes_group.openDataSet("NodesX");
		H5::DataSpace nodes_space = nodes_set.getSpace();

		hsize_t nodecount;
		nodecount = nodes_space.getSimpleExtentNpoints();
		_internalmesh.set_nodecount(nodecount);
		double * nodebuf = new double[nodecount];
		nodes_set.read(nodebuf, H5::PredType::NATIVE_DOUBLE);
		for (auto i=0; i<nodecount; i++) _internalmesh.node(i).set_x(nodebuf[i]);


		HDF5FieldNames.clear();
		H5Literate(nodes_group.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, get_dataset_names_H5, NULL);
		
		for (auto i=0; i<HDF5FieldNames.size(); i++){
			if (HDF5FieldNames.at(i).compare("NodesY") == 0){
				nodes_set = nodes_group.openDataSet("NodesY");
				nodes_set.read(nodebuf, H5::PredType::NATIVE_DOUBLE);
				for (auto i=0; i<nodecount; i++) _internalmesh.node(i).set_y(nodebuf[i]);
			}
			else if(HDF5FieldNames.at(i).compare("NodesZ") == 0){
				nodes_set = nodes_group.openDataSet("NodesZ");
				nodes_set.read(nodebuf, H5::PredType::NATIVE_DOUBLE);
				for (auto i=0; i<nodecount; i++) _internalmesh.node(i).set_z(nodebuf[i]);
			}
		}
		
		_internalmesh.calc_extents();
		delete[] nodebuf;


	// Time reading
	H5::Group time_group = file.openGroup("Time");

		H5::DataSet time_set = time_group.openDataSet("Time");
		H5::DataSpace time_space = time_set.getSpace();
		hsize_t ntime;
		ntime = time_space.getSimpleExtentNpoints();
		_time.resize(ntime);
		
		double timebuf[ntime];
		time_set.read(timebuf, H5::PredType::NATIVE_DOUBLE);
		for (auto i=0; i<ntime; i++) _time[i] = timebuf[i];
		_tstart = timebuf[0]; _tstop = timebuf[ntime-1];
		_dt = timebuf[1] - timebuf[0];
		set_time_span(_tstart, _dt, _tstop);
		//cout << "snapshot size" << _datasnapshots.size() << endl;


	// Fields reading
	H5::Group fields_group = file.openGroup("Fields");

		hsize_t numfields = fields_group.getNumObjs();
		//cout << "there are " << numfields << " fields present in the H5 file" << endl;
		// determine the names of all the datasets in the group
		HDF5FieldNames.clear();
		H5Literate(fields_group.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, get_dataset_names_H5, NULL);
		
		// copy data for each dataset
		//double ** fieldbuf = new double*[nodecount];
		//for (auto i=0; i<nodecount; i++) fieldbuf[i] = new double[ntime];
		double * fieldbuf = new double[nodecount];
		hsize_t fielddims[2];
		fielddims[0] = nodecount;
		fielddims[1] = 1;

		H5::DataSet field_set;
		H5::DataSpace field_space, memspace;
		memspace = H5::DataSpace(2, fielddims);
		hsize_t offset[2], count[2], stride[2], block[2];
		for (auto i=0; i<HDF5FieldNames.size(); i++){
			add_field(HDF5FieldNames.at(i));

			//cout << "working on field: " << HDF5FieldNames.at(i) << endl;
			//cout << "allocated space: " << _datasnapshots.at(0)._datafields.at(HDF5FieldNames.at(i)).size() << endl;

			// open the dataset
			field_set = fields_group.openDataSet(HDF5FieldNames.at(i));
			field_space = field_set.getSpace();

			offset[0] = 0;
			

			count[0]  = nodecount;
			count[1]  = 1;

			stride[0] = 1;
			stride[1] = 1;

			block[0] = 1;
			block[1] = 1;

			//cout << "read field to buffer: " << HDF5FieldNames.at(i) << endl;

			// fill in the data
			for (auto j=0; j<ntime; j++){
				cout << "reading time snapshot " << j << "/" << ntime << "\r" << flush;
				offset[1] = j;
				field_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
				field_set.read(fieldbuf, H5::PredType::NATIVE_DOUBLE, memspace, field_space);

				for (auto k=0; k<nodecount; k++){
					//cout << "j: " << j << " k: " << k << " fieldbuf: " << fieldbuf[k] << endl;// "\r" << flush;
					_datasnapshots.at(j)._datafields.at(HDF5FieldNames.at(i)).at(k) = fieldbuf[k];
				}
			}

		}

		delete[] fieldbuf;
		HDF5FieldNames.clear();

	// close the file
	file.close();
}

// Have to use the C API for this part, because it isn't yet implemented in the C++ API
herr_t get_dataset_names_H5(hid_t loc_id, const char * name, const H5L_info_t *linfo, void *opdata){
	hid_t obj;
	herr_t err;
	H5L_info_t * linkinfo;

	obj = H5Dopen2(loc_id, name, H5P_DEFAULT);

	HDF5FieldNames.push_back(name);
	H5Dclose(obj);
	return 0;
}


SimulationData SimulationData::combine(vector<const SimulationData *> datavec){
	SimulationData combo;
	const Mesh * meshptr;
	double tstart, dt, tstop;
	vector<string> combofields;

	// check that they share the same mesh
	meshptr = datavec.at(0)->_mesh;
	for (auto i=1; i<datavec.size(); i++){
		if (datavec.at(i)->_mesh != meshptr){
			cout << "The input simulations do not share a common mesh!" << endl;
			throw -1;
		}
	}

	// check that the start time is the same
	tstart = datavec.at(0)->_tstart;
	for (auto i=1; i<datavec.size(); i++){
		if (datavec.at(i)->_tstart != tstart){
			cout << "The input simulations do not share a common start time!" << endl;
			throw -1;
		}
	}

	// find the smallest dt and the longest tstop
	dt = datavec.at(0)->_dt;
	tstop = datavec.at(0)->_tstop;
	for (auto i=1; i<datavec.size(); i++){
		if (datavec.at(i)->_dt < dt) dt = datavec.at(i)->_dt;
		if (datavec.at(i)->_tstop > tstop) tstop = datavec.at(i)->_tstop;
	}


	// checking is done... bind stuff
	combo.bind_mesh(*meshptr);
	combo.set_time_span(tstart, dt, tstop);

	// find unique fieldnames
	for (auto i=0; i<datavec.size(); i++){
		combofields.insert(combofields.end(), datavec.at(i)->_fieldnames.begin(), datavec.at(i)->_fieldnames.end());
	}
	sort(combofields.begin(), combofields.end());
	combofields.erase(unique(combofields.begin(), combofields.end()), combofields.end());

	for (auto i=0; i<combofields.size(); i++) combo.add_field(combofields.at(i));

	// populate the fields for each fieldname
	double tcur = tstart;
	unsigned int tcur_ind = 0;
	for (auto i=0; i<combofields.size(); i++){
		cout << "on field: " << combofields.at(i) << endl;
		for (auto j=0; j<datavec.size(); j++){
			if (!datavec.at(j)->field_present(combofields.at(i))) continue;

			cout << "field is present for datavec #" << j << endl;
			// if the field is present, copy it for each time slot
			// some could have different time steps though, so check for that
			for (auto k=0; k<combo._time.size(); k++){
				combo.add_data_at_index(k, combofields.at(i), datavec.at(j)->get_data_at_time(tcur, combofields.at(i)));
				if (datavec.at(j)->_dt == dt) tcur_ind++;
				else if (tcur + datavec.at(j)->_dt >= tstart + (k+1)*dt) cout << "inc+" << tcur_ind++; 
				if (tcur_ind > datavec.at(j)->_time.size()-1) tcur_ind--;
				cout << "filled in k: " << k << "/" << combo._time.size() << "\r" << flush;
				tcur = datavec.at(j)->_time.at(tcur_ind);
			}
			cout << endl;
			cout << "finished copying field" << endl;
			cout << endl;

			// once the field is copied, move on to the next one
			break;
		}
	}

	combo.print_summary();
	return combo;
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

