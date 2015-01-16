// this is a placeholder for the PointCloud class
//#define _TEST_
// it should hold the following elements (at least):
//  - x,y,z,intensity
//  - return number
//  - classification
//  - extra fields, extra field names
//  - projection (UTM, ECEF, lat/lon) (choose one to always work in... I'm thinking ECEF)
//  - pointcount

// and it should have the following functions (at least):
//  - read LAS file (DONE)
//  - write LAS file
//  - read ascii delimited x,y,z,i,c
//  - write ascii delimited x,y,z,i,c
//  - print a summary of data (DONE)
//  - recalculate the x,y,z extents (DONE)
//  - conversion to/from ECEF/latlon/UTM
//  - create/destroy data vectors (DONE)
#include "PointCloud.hpp"

using namespace std;

PointCloud::PointCloud(unsigned int numpts){
  pointcount = numpts;
  x = new double[pointcount];
  y = new double[pointcount];
  z = new double[pointcount];
  gpstime = NULL;
  intensity = NULL;
  classification = NULL;
  RGB = NULL;

  xmin = 0; xmax = 0; ymin = 0; ymax = 0; zmin = 0; zmax = 0;
  gpst_min = 0; gpst_max = 0;
}

PointCloud::~PointCloud(){
  delete[] x;
  delete[] y;
  delete[] z;
  if (gpstime != NULL) delete[] gpstime;
  if (intensity != NULL) delete[] intensity;
  if (classification != NULL) delete[] classification;
  if (RGB != NULL) delete[] RGB;
}

void PointCloud::print_summary(){
  cout << "Point Cloud Summary:" << endl;
  cout << "  pointcount: " << pointcount << endl;
  cout << "  existing fields: x" << endl;
  cout << "                   y" << endl;
  cout << "                   z" << endl;
  if (gpstime !=NULL) cout << "                   gpstime" << endl;
  if (intensity !=NULL) cout << "                   intensity" << endl;
  if (classification !=NULL) cout << "                   classification" << endl;
  if (RGB !=NULL) cout << "                   RGB" << endl;
  cout << "  data extents:" << endl;
  cout << "       x:[" << xmin << ", " << xmax << "]" << endl;
  cout << "       y:[" << ymin << ", " << ymax << "]" << endl;
  cout << "       z:[" << zmin << ", " << zmax << "]" << endl;
  if (gpstime != NULL) cout << "       gpstime:[" << gpst_min << ", " << gpst_max << "]" << endl;

  return;
}

void PointCloud::calc_extents(){
  xmin = x[0]; xmax = x[0]; ymin = y[0]; ymax = y[0];
  zmin = z[0]; zmax = z[0];
  for (unsigned int i=1; i<pointcount; i++){
    if (x[i] < xmin) xmin = x[i];
    else if (x[i] > xmax) xmax = x[i];
    if (y[i] < ymin) ymin = y[i];
    else if (y[i] > ymax) ymax = y[i];
    if (z[i] < zmin) zmin = z[i];
    else if (z[i] > zmax) zmax = z[i];
  }
  if (gpstime != NULL){
    gpst_min = gpstime[0]; gpst_max = gpstime[0];
    for (unsigned int i=1; i<pointcount; i++){
      if (gpstime[i] < gpst_min) gpst_min = gpstime[i];
      else if (gpstime[i] > gpst_max) gpst_max = gpstime[i];
    }
  }

  return;
}

void PointCloud::add_intensity(){
  if (intensity != NULL) cout << "intensity already exists!" << endl;
  else intensity = new unsigned short[pointcount];
}

void PointCloud::add_classification(){
  if (classification != NULL) cout << "classification already exists!" << endl;
  else classification = new unsigned char[pointcount];
}

void PointCloud::add_gpstime(){
  if (gpstime != NULL) cout << "gpstime already exists!" << endl;
  else gpstime = new double[pointcount];
}

void PointCloud::add_RGB(){
  if (RGB != NULL) cout << "RGB already exists!" << endl;
  else RGB = new rgb48[pointcount];
}


PointCloud * PointCloud::read_LAS(string filename, unsigned int byte_offset){
  // define vars
  bool fieldexist=false;
  char signature[4];
  unsigned char ver_major, ver_minor, point_format_id, dummy_1;
  unsigned short header_size, point_record_bytes, dummy_2;
  unsigned int point_offset, dummy_4;
  unsigned int pt_count;
  double x_scale, y_scale, z_scale, x_offset, y_offset, z_offset;
  double x_min, y_min, z_min, x_max, y_max, z_max, dummy_8;
  int fd;
  PointCloud * cloud;
  char * lasmap;
  size_t res;

  // open the file
  FILE * fid;
  fid = fopen(filename.c_str(), "rb");

  // do a size check
  fseek(fid, 0L, SEEK_END);
  size_t sz = ftell(fid);
  fseek(fid, 0L, SEEK_SET);
  //cout << "the file is " << sz << " bytes long" << endl;
  
  // move to byte_offset
  fseek(fid, byte_offset, SEEK_SET);

  // check LASF signature
  res = fread(&signature, 4, 1, fid); // "LASF" signature

  //cout << "signature: " << signature[0] << signature[1] << signature[2] << signature[3] << endl;

  // collect header info
  res = fread(&dummy_2, 2, 1, fid); 		// File Source ID
  res = fread(&dummy_2, 2, 1, fid); 		// Global Encoding
  res = fread(&dummy_4, 4, 1, fid); 		// Project ID - GUID data 1
  res = fread(&dummy_2, 2, 1, fid); 		// Project ID - GUID data 2
  res = fread(&dummy_2, 2, 1, fid); 		// Project ID - GUID data 3
  res = fread(&dummy_8, 8, 1, fid); 		// Project ID - GUID data 4
  res = fread(&ver_major, 1, 1, fid); 	// Version Major
  res = fread(&ver_minor, 1, 1, fid);   // Version Minor
  for (int i=0; i<32; i++) res = fread(&dummy_1, 1, 1, fid); // System Identifier
  for (int i=0; i<32; i++) res = fread(&dummy_1, 1, 1, fid); // Generating Software
  res = fread(&dummy_2, 2, 1, fid); 		// File Creation Day of Year
  res = fread(&dummy_2, 2, 1, fid); 		// File Creation Year
  res = fread(&header_size, 2, 1, fid); // Header Size
  res = fread(&point_offset, 4, 1, fid);// Offset to Point Data
  res = fread(&dummy_4, 4, 1, fid); 		// Number of Variable Length Records
  res = fread(&point_format_id, 1, 1, fid);	// Point Data Format ID
  res = fread(&point_record_bytes, 2, 1, fid);	// Point Data Record Length
  res = fread(&pt_count, 4, 1, fid);		// Number of point records
  for (int i=0; i<5; i++) res = fread(&dummy_4, 4, 1, fid); // Number of points by return
  if (ver_minor > 2) res = fread(&dummy_8, 8, 1, fid); // extra points by return 
  res = fread(&x_scale, 8, 1, fid);		  // X Scale Factor
  res = fread(&y_scale, 8, 1, fid);		  // Y Scale Factor
  res = fread(&z_scale, 8, 1, fid);		  // Z Scale Factor
  res = fread(&x_offset, 8, 1, fid);		// X Offset
  res = fread(&y_offset, 8, 1, fid);		// Y Offset
  res = fread(&z_offset, 8, 1, fid);		// Z Offset
  res = fread(&x_max, 8, 1, fid);		    // X Max
  res = fread(&x_min, 8, 1, fid);		    // X Min
  res = fread(&y_max, 8, 1, fid);		    // Y Max
  res = fread(&y_min, 8, 1, fid);		    // Y Min
  res = fread(&z_max, 8, 1, fid);		    // Z Max
  res = fread(&z_min, 8, 1, fid);		    // Z Min

  // verify that the number of points match up
  //cout << "we have mapped" << point_offset + point_record_bytes*pt_count << " bytes" << endl;
  //cout << "the file should have " << (sz-point_offset)/point_record_bytes << " points" << endl;
  unsigned int realsize = (sz-point_offset - byte_offset)/point_record_bytes;
  if (pt_count > realsize){
    cout << "WARNING: byte count doesn't match up with reported point count" << endl;
    cout << "WARNING: proceeding with calculated byte count" << endl;
    pt_count = realsize;
  }

  // read projection info from Variable Length Records (VLRs)

  // close the file
  fclose(fid);

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

  /*
  cout << "LAS version " << int(ver_major) << "." << int(ver_minor) << endl;
  cout << "point format id: " << int(point_format_id) << endl;
  cout << "bytes in point record: " << point_record_bytes << endl;
  cout << "pointcount: " << pt_count << endl;
  cout << "point offset: " << point_offset << endl;
  cout << "xoffset: " << x_offset << " yoffset:" << y_offset << " zoffset:" << z_offset << endl;
  cout << "xscale: " << x_scale << " yscale:" << y_scale << " zscale:" << z_scale << endl;
  cout << "xmax: " << x_max << " xmin: " << x_min << endl;
  cout << "ymax: " << y_max << " ymin: " << y_min << endl;
  cout << "zmax: " << z_max << " zmin: " << z_min << endl;
  */

// initialize point cloud
  cloud = new PointCloud(pt_count);
  cloud->add_intensity();
  cloud->add_classification();
  cloud->xmin = x_min;
  cloud->xmax = x_max;
  cloud->ymin = y_min;
  cloud->ymax = y_max;
  cloud->zmin = z_min;
  cloud->zmax = z_max;

  cloud->read_LAS_internal(filename, byte_offset);

  return cloud;

/*
  // memory map the point data
  fd = open(filename, O_RDONLY);
  if (fd < 0){
    cout << "Error opening file in readLas" << endl;
    throw -1;
  }
  lseek(fd, byte_offset, SEEK_SET);
  lasmap = (char *)mmap(NULL, point_offset + point_record_bytes*pt_count, PROT_READ, MAP_PRIVATE, fd, 0);
  if (lasmap == MAP_FAILED){
    cout << "Failed to map las file" << endl;
    throw -1;
  }  


  // do a switch for point record format and memcpy the points
  switch (point_format_id)
  {
    case 0:
      // initialize array of structs
      las_pt_0 *laspts0;
      laspts0 = new las_pt_0[pt_count];

      //cout << "about to memcpy " << sizeof(las_pt_0)*pt_count << " " << " " << point_record_bytes*pt_count << endl;

      // copy from memory map into array of structs
      memcpy(laspts0, &lasmap[point_offset], point_record_bytes*pt_count);

      //cout << "about to extract points" << endl;
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts0[i].X)*x_scale + x_offset; 		// X
        cloud->y[i] = double(laspts0[i].Y)*y_scale + y_offset;		// Y
        cloud->z[i] = double(laspts0[i].Z)*z_scale + z_offset;		// Z
        cloud->intensity[i] = laspts0[i].Intensity;			// Intensity
        cloud->classification[i] = laspts0[i].Classification;		// Classification
      }

      delete[] laspts0;
      //cout << "successfully deleted" << endl;
      break;

    case 1:
      // initialize array of structs
      las_pt_1 *laspts1;
      laspts1 = new las_pt_1[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts1, &lasmap[point_offset-1], point_record_bytes*pt_count);

      cloud->add_gpstime();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts1[i].X)*x_scale + x_offset;    // X
        cloud->y[i] = double(laspts1[i].Y)*y_scale + y_offset;    // Y
        cloud->z[i] = double(laspts1[i].Z)*z_scale + z_offset;    // Z
        cloud->intensity[i] = laspts1[i].Intensity;                // Intensity
        cloud->classification[i] = laspts1[i].Classification;      // Classification
        cloud->gpstime[i] = laspts1[i].GPSTime;
      }
      cloud->gpst_min = cloud->gpstime[0];
      cloud->gpst_max = cloud->gpstime[pt_count-1];

      delete[] laspts1;
      break;

    case 2:
      // initialize array of structs
      las_pt_2 *laspts2;
      laspts2 = new las_pt_2[pt_count];

      //cout << "about to memcpy" << sizeof(las_pt_2)*pt_count << " " << " " << point_record_bytes*pt_count << endl;
      // copy from memory map into array of structs
      memcpy(laspts2, &lasmap[point_offset], point_record_bytes*pt_count);

      //cout << "about to copy stuff" << endl;
      cloud->add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts2[i].X)*x_scale + x_offset;    // X
        cloud->y[i] = double(laspts2[i].Y)*y_scale + y_offset;    // Y
        cloud->z[i] = double(laspts2[i].Z)*z_scale + z_offset;    // Z
        cloud->intensity[i] = laspts2[i].Intensity;                // Intensity
        cloud->classification[i] = laspts2[i].Classification;      // Classification
        cloud->RGB[i].R = laspts2[i].Red;
        cloud->RGB[i].G = laspts2[i].Green;
        cloud->RGB[i].B = laspts2[i].Blue;
      }

      delete[] laspts2;
      break;

    case 3:
      // initialize array of structs
      las_pt_3 *laspts3;
      laspts3 = new las_pt_3[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts3, &lasmap[point_offset-1], point_record_bytes*pt_count);

      cloud->add_gpstime();
      cloud->add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts3[i].X)*x_scale + x_offset;    // X
        cloud->y[i] = double(laspts3[i].Y)*y_scale + y_offset;    // Y
        cloud->z[i] = double(laspts3[i].Z)*z_scale + z_offset;    // Z
        cloud->intensity[i] = laspts3[i].Intensity;                // Intensity
        cloud->classification[i] = laspts3[i].Classification;      // Classification
        cloud->gpstime[i] = laspts3[i].GPSTime;
        cloud->RGB[i].R = laspts3[i].Red;
        cloud->RGB[i].G = laspts3[i].Green;
        cloud->RGB[i].B = laspts3[i].Blue;
      }
      cloud->gpst_min = cloud->gpstime[0];
      cloud->gpst_max = cloud->gpstime[pt_count-1];

      delete[] laspts3;
      break;

    case 4:
      // initialize array of structs
      las_pt_4 *laspts4;
      laspts4 = new las_pt_4[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts4, &lasmap[point_offset-1], point_record_bytes*pt_count);

      cloud->add_gpstime();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts4[i].X)*x_scale + x_offset;    // X
        cloud->y[i] = double(laspts4[i].Y)*y_scale + y_offset;    // Y
        cloud->z[i] = double(laspts4[i].Z)*z_scale + z_offset;    // Z
        cloud->intensity[i] = laspts4[i].Intensity;                // Intensity
        cloud->classification[i] = laspts4[i].Classification;      // Classification
        cloud->gpstime[i] = laspts4[i].GPSTime;
      }
      cloud->gpst_min = cloud->gpstime[0];
      cloud->gpst_max = cloud->gpstime[pt_count-1];

      delete[] laspts4;
      break;

    case 5:
      // initialize array of structs
      las_pt_5 *laspts5;
      laspts5 = new las_pt_5[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts5, &lasmap[point_offset-1], point_record_bytes*pt_count);

      cloud->add_gpstime();
      cloud->add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        cloud->x[i] = double(laspts5[i].X)*x_scale + x_offset;    // X
        cloud->y[i] = double(laspts5[i].Y)*y_scale + y_offset;    // Y
        cloud->z[i] = double(laspts5[i].Z)*z_scale + z_offset;    // Z
        cloud->intensity[i] = laspts5[i].Intensity;                // Intensity
        cloud->classification[i] = laspts5[i].Classification;      // Classification
        cloud->gpstime[i] = laspts5[i].GPSTime;
        cloud->RGB[i].R = laspts5[i].Red;
        cloud->RGB[i].G = laspts5[i].Green;
        cloud->RGB[i].B = laspts5[i].Blue;
      }
      cloud->gpst_min = cloud->gpstime[0];
      cloud->gpst_max = cloud->gpstime[pt_count-1];

      delete[] laspts5;
      break;

    default:
      cout << "this should not be possible" << endl;
  }

  // unmap the data
    //cout << "trying to unmap" << endl;

  if (munmap(lasmap, point_offset + point_record_bytes*pt_count) < 0){
    cout << "ruh roh! problem unmapping LAS file" << endl;
    throw -1;
  }
  close(fd);
  //cout << "finished unmapping" << endl;
  

  // check to see intensity and classification contain actual info
  for (unsigned int i=0; i<pt_count; i++){
    if (cloud->intensity[i] != 0) {
      fieldexist = true;
      break;
    }
  }
  if (!fieldexist){
    delete[] cloud->intensity;
    cloud->intensity = NULL;
  }
  fieldexist = false;
  for (unsigned int i=0; i<pt_count; i++){
    if (cloud->classification[i] != 0) {
      fieldexist = true;
      break;
    }
  }
  if (!fieldexist){
    delete[] cloud->classification;
    cloud->classification = NULL;
  }
  fieldexist = false;

  return cloud;
}
*/
}

void PointCloud::read_LAS_internal(string filename, unsigned int byte_offset){
  // define vars
  bool fieldexist=false;
  char signature[4];
  unsigned char ver_major, ver_minor, point_format_id, dummy_1;
  unsigned short header_size, point_record_bytes, dummy_2;
  unsigned int point_offset, dummy_4;
  unsigned int pt_count;
  double x_scale, y_scale, z_scale, x_offset, y_offset, z_offset;
  double x_min, y_min, z_min, x_max, y_max, z_max, dummy_8;
  int fd;
  char * lasmap;
  size_t res;

  // open the file
  FILE * fid;
  fid = fopen(filename.c_str(), "rb");

  // do a size check
  fseek(fid, 0L, SEEK_END);
  size_t sz = ftell(fid);
  fseek(fid, 0L, SEEK_SET);
  //cout << "the file is " << sz << " bytes long" << endl;
  
  // move to byte_offset
  fseek(fid, byte_offset, SEEK_SET);

  // check LASF signature
  res = fread(&signature, 4, 1, fid); // "LASF" signature

  //cout << "signature: " << signature[0] << signature[1] << signature[2] << signature[3] << endl;

  // collect header info
  res = fread(&dummy_2, 2, 1, fid);     // File Source ID
  res = fread(&dummy_2, 2, 1, fid);     // Global Encoding
  res = fread(&dummy_4, 4, 1, fid);     // Project ID - GUID data 1
  res = fread(&dummy_2, 2, 1, fid);     // Project ID - GUID data 2
  res = fread(&dummy_2, 2, 1, fid);     // Project ID - GUID data 3
  res = fread(&dummy_8, 8, 1, fid);     // Project ID - GUID data 4
  res = fread(&ver_major, 1, 1, fid);   // Version Major
  res = fread(&ver_minor, 1, 1, fid);   // Version Minor
  for (int i=0; i<32; i++) res = fread(&dummy_1, 1, 1, fid); // System Identifier
  for (int i=0; i<32; i++) res = fread(&dummy_1, 1, 1, fid); // Generating Software
  res = fread(&dummy_2, 2, 1, fid);     // File Creation Day of Year
  res = fread(&dummy_2, 2, 1, fid);     // File Creation Year
  res = fread(&header_size, 2, 1, fid); // Header Size
  res = fread(&point_offset, 4, 1, fid);// Offset to Point Data
  res = fread(&dummy_4, 4, 1, fid);     // Number of Variable Length Records
  res = fread(&point_format_id, 1, 1, fid); // Point Data Format ID
  res = fread(&point_record_bytes, 2, 1, fid);  // Point Data Record Length
  res = fread(&pt_count, 4, 1, fid);    // Number of point records
  for (int i=0; i<5; i++) res = fread(&dummy_4, 4, 1, fid); // Number of points by return
  if (ver_minor > 2) res = fread(&dummy_8, 8, 1, fid); // extra points by return 
  res = fread(&x_scale, 8, 1, fid);     // X Scale Factor
  res = fread(&y_scale, 8, 1, fid);     // Y Scale Factor
  res = fread(&z_scale, 8, 1, fid);     // Z Scale Factor
  res = fread(&x_offset, 8, 1, fid);    // X Offset
  res = fread(&y_offset, 8, 1, fid);    // Y Offset
  res = fread(&z_offset, 8, 1, fid);    // Z Offset
  res = fread(&x_max, 8, 1, fid);       // X Max
  res = fread(&x_min, 8, 1, fid);       // X Min
  res = fread(&y_max, 8, 1, fid);       // Y Max
  res = fread(&y_min, 8, 1, fid);       // Y Min
  res = fread(&z_max, 8, 1, fid);       // Z Max
  res = fread(&z_min, 8, 1, fid);       // Z Min

  // verify that the number of points match up
  //cout << "we have mapped" << point_offset + point_record_bytes*pt_count << " bytes" << endl;
  //cout << "the file should have " << (sz-point_offset)/point_record_bytes << " points" << endl;
  unsigned int realsize = (sz-point_offset - byte_offset)/point_record_bytes;
  if (pt_count > realsize){
    cout << "WARNING: byte count doesn't match up with reported point count" << endl;
    cout << "WARNING: proceeding with calculated byte count" << endl;
    pt_count = realsize;
  }

  // read projection info from Variable Length Records (VLRs)

  // close the file
  fclose(fid);

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

  /*
  cout << "LAS version " << int(ver_major) << "." << int(ver_minor) << endl;
  cout << "point format id: " << int(point_format_id) << endl;
  cout << "bytes in point record: " << point_record_bytes << endl;
  cout << "pointcount: " << pt_count << endl;
  cout << "point offset: " << point_offset << endl;
  cout << "xoffset: " << x_offset << " yoffset:" << y_offset << " zoffset:" << z_offset << endl;
  cout << "xscale: " << x_scale << " yscale:" << y_scale << " zscale:" << z_scale << endl;
  cout << "xmax: " << x_max << " xmin: " << x_min << endl;
  cout << "ymax: " << y_max << " ymin: " << y_min << endl;
  cout << "zmax: " << z_max << " zmin: " << z_min << endl;
  */

  // memory map the point data
  fd = open(filename.c_str(), O_RDONLY);
  if (fd < 0){
    cout << "Error opening file in readLas" << endl;
    throw -1;
  }
  lseek(fd, byte_offset, SEEK_SET);
  lasmap = (char *)mmap(NULL, point_offset + point_record_bytes*pt_count, PROT_READ, MAP_PRIVATE, fd, 0);
  if (lasmap == MAP_FAILED){
    cout << "Failed to map las file" << endl;
    throw -1;
  }  


  // do a switch for point record format and memcpy the points
  switch (point_format_id)
  {
    case 0:
      // initialize array of structs
      las_pt_0 *laspts0;
      laspts0 = new las_pt_0[pt_count];

      //cout << "about to memcpy " << sizeof(las_pt_0)*pt_count << " " << " " << point_record_bytes*pt_count << endl;

      // copy from memory map into array of structs
      memcpy(laspts0, &lasmap[point_offset], point_record_bytes*pt_count);

      //cout << "about to extract points" << endl;
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts0[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts0[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts0[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts0[i].Intensity;     // Intensity
        classification[i] = laspts0[i].Classification;   // Classification
      }

      delete[] laspts0;
      //cout << "successfully deleted" << endl;
      break;

    case 1:
      // initialize array of structs
      las_pt_1 *laspts1;
      laspts1 = new las_pt_1[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts1, &lasmap[point_offset-1], point_record_bytes*pt_count);

      add_gpstime();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts1[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts1[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts1[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts1[i].Intensity;                // Intensity
        classification[i] = laspts1[i].Classification;      // Classification
        gpstime[i] = laspts1[i].GPSTime;
      }
      gpst_min = gpstime[0];
      gpst_max = gpstime[pt_count-1];

      delete[] laspts1;
      break;

    case 2:
      // initialize array of structs
      las_pt_2 *laspts2;
      laspts2 = new las_pt_2[pt_count];

      //cout << "about to memcpy" << sizeof(las_pt_2)*pt_count << " " << " " << point_record_bytes*pt_count << endl;
      // copy from memory map into array of structs
      memcpy(laspts2, &lasmap[point_offset], point_record_bytes*pt_count);

      //cout << "about to copy stuff" << endl;
      add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts2[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts2[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts2[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts2[i].Intensity;                // Intensity
        classification[i] = laspts2[i].Classification;      // Classification
        RGB[i].R = laspts2[i].Red;
        RGB[i].G = laspts2[i].Green;
        RGB[i].B = laspts2[i].Blue;
      }

      delete[] laspts2;
      break;

    case 3:
      // initialize array of structs
      las_pt_3 *laspts3;
      laspts3 = new las_pt_3[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts3, &lasmap[point_offset-1], point_record_bytes*pt_count);

      add_gpstime();
      add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts3[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts3[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts3[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts3[i].Intensity;                // Intensity
        classification[i] = laspts3[i].Classification;      // Classification
        gpstime[i] = laspts3[i].GPSTime;
        RGB[i].R = laspts3[i].Red;
        RGB[i].G = laspts3[i].Green;
        RGB[i].B = laspts3[i].Blue;
      }
      gpst_min = gpstime[0];
      gpst_max = gpstime[pt_count-1];

      delete[] laspts3;
      break;

    case 4:
      // initialize array of structs
      las_pt_4 *laspts4;
      laspts4 = new las_pt_4[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts4, &lasmap[point_offset-1], point_record_bytes*pt_count);

      add_gpstime();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts4[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts4[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts4[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts4[i].Intensity;                // Intensity
        classification[i] = laspts4[i].Classification;      // Classification
        gpstime[i] = laspts4[i].GPSTime;
      }
      gpst_min = gpstime[0];
      gpst_max = gpstime[pt_count-1];

      delete[] laspts4;
      break;

    case 5:
      // initialize array of structs
      las_pt_5 *laspts5;
      laspts5 = new las_pt_5[pt_count];

      // copy from memory map into array of structs
      memcpy(laspts5, &lasmap[point_offset-1], point_record_bytes*pt_count);

      add_gpstime();
      add_RGB();
      // loop over structs to extract the points
      for (unsigned int i=0; i<pt_count; i++){
        x[i] = double(laspts5[i].X)*x_scale + x_offset;    // X
        y[i] = double(laspts5[i].Y)*y_scale + y_offset;    // Y
        z[i] = double(laspts5[i].Z)*z_scale + z_offset;    // Z
        intensity[i] = laspts5[i].Intensity;                // Intensity
        classification[i] = laspts5[i].Classification;      // Classification
        gpstime[i] = laspts5[i].GPSTime;
        RGB[i].R = laspts5[i].Red;
        RGB[i].G = laspts5[i].Green;
        RGB[i].B = laspts5[i].Blue;
      }
      gpst_min = gpstime[0];
      gpst_max = gpstime[pt_count-1];

      delete[] laspts5;
      break;

    default:
      cout << "this should not be possible" << endl;
  }

  // unmap the data
    //cout << "trying to unmap" << endl;

  if (munmap(lasmap, point_offset + point_record_bytes*pt_count) < 0){
    cout << "ruh roh! problem unmapping LAS file" << endl;
    throw -1;
  }
  close(fd);
  //cout << "finished unmapping" << endl;
  

  // check to see intensity and classification contain actual info
  for (unsigned int i=0; i<pt_count; i++){
    if (intensity[i] != 0) {
      fieldexist = true;
      break;
    }
  }
  if (!fieldexist){
    delete[] intensity;
    intensity = NULL;
  }
  fieldexist = false;
  for (unsigned int i=0; i<pt_count; i++){
    if (classification[i] != 0) {
      fieldexist = true;
      break;
    }
  }
  if (!fieldexist){
    delete[] classification;
    classification = NULL;
  }
  fieldexist = false;

  return;
}


void PointCloud::write_LAS(string filename){
  cout << "writing not quite supported yet" << endl;
}

PointCloud * PointCloud::subset(bool * keep){
  // declare vars
  unsigned int subset_count=0, ct;
  PointCloud * cloud_subset;

  for (unsigned int i=0; i<pointcount; i++){
    if (keep[i] == true) subset_count++;
  }

  cloud_subset = new PointCloud(subset_count);

  // copy x, y, z data
  ct=0;
  for (unsigned int i=0; i<pointcount; i++){
    if (keep[i]){
      cloud_subset->x[ct] = x[i];
      cloud_subset->y[ct] = y[i];
      cloud_subset->z[ct] = z[i];
      ct++;
    }
  }

  // copy intensity data
  ct=0;
  if (intensity != NULL){
    cloud_subset->add_intensity();
    for (unsigned int i=0; i<pointcount; i++){
      if (keep[i]){
        cloud_subset->intensity[ct] = intensity[i];
        ct++;
      }
    }
  }

  // copy classification data
  ct=0;
  if (classification != NULL){
    cloud_subset->add_classification();
    for (unsigned int i=0; i<pointcount; i++){
      if (keep[i]){
        cloud_subset->classification[ct] = classification[i];
        ct++;
      }
    }
  }

  // copy gpstime data
  ct=0;
  if (gpstime != NULL){
    cloud_subset->add_gpstime();
    for (unsigned int i=0; i<pointcount; i++){
      if (keep[i]){
        cloud_subset->gpstime[ct] = gpstime[i];
        ct++;
      }
    }
  }

  // copy RGB data
  ct=0;
  if (RGB != NULL){
    cloud_subset->add_RGB();
    for (unsigned int i=0; i<pointcount; i++){
      if (keep[i]){
        cloud_subset->RGB[ct].R = RGB[i].R;
        cloud_subset->RGB[ct].G = RGB[i].G;
        cloud_subset->RGB[ct].B = RGB[i].B;
        ct++;
      }
    }
  }

  cloud_subset->calc_extents();

  return cloud_subset;
}

PointCloud * PointCloud::subset(unsigned int * keep_inds, unsigned int keep_count){
  // declare vars
  bool * keep;
  PointCloud * cloud_subset;

  keep = new bool[pointcount];

  for (unsigned int i=0; i<pointcount; i++) keep[i] = false;
  for (unsigned int j=0; j<keep_count; j++) keep[keep_inds[j]] = true;

  cloud_subset = subset(keep);

  delete[] keep;

  return cloud_subset;
}


#ifdef _TEST_

int main(int argc, char * argv[]){
  // declare vars
  PointCloud * cloud, * cloud_sub;
  unsigned int * keep_inds;

  // take one command line argument as the las file name and read it
  cloud = PointCloud::read_LAS(argv[1]);

  // output summary info
  cloud->print_summary();

  cloud->calc_extents();
  cloud->print_summary();

  // test the subset function
  keep_inds = new unsigned int[4];
  keep_inds[0] = 1;
  keep_inds[1] = 100;
  keep_inds[2] = 200;
  keep_inds[3] = 300;
  cloud_sub = cloud->subset(keep_inds, 4);
  cloud_sub->print_summary();

  // delete stuff
  delete cloud;
  delete cloud_sub;

  return 0;
}
#endif