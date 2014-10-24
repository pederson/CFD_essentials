// this is a placeholder for the PointCloud class

// it should hold the following elements (at least):
//  - x,y,z,intensity
//  - return number
//  - classification
//  - extra fields, extra field names
//  - projection (UTM, ECEF, lat/lon) (choose one to always work in... I'm thinking ECEF)
//  - pointcount

// and it should have the following functions (at least):
//  - read LAS file
//  - write LAS file
//  - read ascii delimited x,y,z,i,c
//  - write ascii delimited x,y,z,i,c
//  - print a summary of data
//  - recalculate the x,y,z extents
//  - conversion to/from ECEF/latlon/UTM
//  - create/destroy data vectors
#include "PointCloud.hpp"

using namespace std;

PointCloud * PointCloud::readLas(char * filename, unsigned int byte_offset){
  // define vars
  char signature[4];
  unsigned char ver_major, ver_minor, point_format_id, dummy_1;
  unsigned short header_size, point_record_bytes, dummy_2;
  unsigned long point_offset, dummy_4;
  unsigned int pt_count;
  double x_scale, y_scale, z_scale, x_offset, y_offset, z_offset;
  double x_min, y_min, z_min, x_max, y_max, z_max, dummy_8;
  int fd;

  // open the file
  ifstream lasfile;
  lasfile.open(filename);
  
  // move to byte_offset
  lasfile.seekg(byte_offset);

  // check LASF signature
  lasfile.read(signature, 4); // "LASF" signature

  // collect header info
  lasfile.read(&dummy_2, 2); 		// File Source ID
  lasfile.read(&dummy_2, 2); 		// Global Encoding
  lasfile.read(&dummy_4, 4); 		// Project ID - GUID data 1
  lasfile.read(&dummy_2, 2); 		// Project ID - GUID data 2
  lasfile.read(&dummy_2, 2); 		// Project ID - GUID data 3
  lasfile.read(&dummy_8, 8); 		// Project ID - GUID data 4
  lasfile.read(&ver_major, 1); 		// Version Major
  lasfile.read(&ver_minor, 1); 		// Version Minor
  for (int i=0; i<32; i++) lasfile.read(&dummy_1, 1); // System Identifier
  for (int i=0; i<32; i++) lasfile.read(&dummy_1, 1); // Generating Software
  lasfile.read(&dummy_2, 2); 		// File Creation Day of Year
  lasfile.read(&dummy_2, 2); 		// File Creation Year
  lasfile.read(&header_size, 2); 	// Header Size
  lasfile.read(&point_offset, 4); 	// Offset to Point Data
  lasfile.read(&dummy_4, 4); 		// Number of Variable Length Records
  lasfile.read(&point_format_id, 1);	// Point Data Format ID
  lasfile.read(&point_record_bytes, 2);	// Point Data Record Length
  lasfile.read(&pt_count, 4);		// Number of point records
  for (int i=0; i<7; i++) lasfile.read(&dummy_4, 4); // Number of points by return
  lasfile.read(&x_scale, 8);		// X Scale Factor
  lasfile.read(&y_scale, 8);		// Y Scale Factor
  lasfile.read(&z_scale, 8);		// Z Scale Factor
  lasfile.read(&x_offset, 8);		// X Offset
  lasfile.read(&y_offset, 8);		// Y Offset
  lasfile.read(&z_offset, 8);		// Z Offset
  lasfile.read(&x_max, 8);		// X Max
  lasfile.read(&x_min, 8);		// X Min
  lasfile.read(&y_max, 8);		// Y Max
  lasfile.read(&y_min, 8);		// Y Min
  lasfile.read(&z_max, 8);		// Z Max
  lasfile.read(&z_min, 8);		// Z Min

  // close the file
  lasfile.close();

  // check the version number and alert the user if it is not fully supported
  if (ver_minor > 2 || ver_major > 1){
    cout << "Some features of LAS version " << ver_major << "." << ver_minor ;
    cout << " may not be supported" << endl;
  }

  // check the record format number and alert the user if it is not supported
  if (point_format_id > 3){
    cout << "Some features of point format " << point_format_id ;
    cout << " may not be supported" << endl;
  }
  if (point_format_id > 5){
    cout << "Point format " << point_format_id << " is not supported" << endl;
    throw -1;
  }

  // memory map the point data
  fd = open(filename, O_RDONLY);
  if (fd == -1){
    cout << "Error opening file in readLas" << endl;
    throw -1;
  }
  lasmap = mmap(NULL, point_record_bytes*pt_count + point_offset, PROT_READ, MAP_PRIVATE, fd, point_offset);
  if (lasmap == MAP_FAILED){
    cout << "Failed to map las file" << endl;
    throw -1;
  }
  close(fd);
  


  // do a switch for point record format and memcpy the points


  // unmap the data


  // delete stuff

  return cloud;
}

void PointCloud::writeLas(char * filename){

}
