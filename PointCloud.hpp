#ifndef _POINTCLOUD_H
#define _POINTCLOUD_H

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include <sys/mman.h>

// this is a placeholder for the PointCloud class

// it should hold the following elements (at least):
//  - x,y,z,intensity
//  - return number
//  - classification
//  - extra fields, extra field names
//  - projection (UTM, ECEF, lat/lon) (choose one to always work in... I'm thinking ECEF)
//
// and it should have the following functions (at least):
//  - read LAS file
//  - write LAS file
//  - read ascii delimited x,y,z,i,c
//  - write ascii delimited x,y,z,i,c
//  - print a summary of data
//  - recalculate the x,y,z extents
//  - conversion to/from ECEF/latlon/UTM
//  - create/destroy data vectors


#pragma pack(push,1)
struct rgb48{
  unsigned short R;
  unsigned short G;
  unsigned short B;
};
#pragma pack(pop)

class PointCloud{
public:

  PointCloud(unsigned int numpts);
  ~PointCloud();

  // base member data accessors
  unsigned int size() const {return pointcount;};
  double x_max() const {return xmax;};
  double x_min() const {return xmin;};
  double y_max() const {return ymax;};
  double y_min() const {return ymin;};
  double z_max() const {return zmax;};
  double z_min() const {return zmin;};
  const double * x_ptr() const {return x;};
  const double * y_ptr() const {return y;};
  const double * z_ptr() const {return z;};
  //double x(unsigned int i) {return x[i];};
  //double y(unsigned int i) {return y[i];};
  //double z(unsigned int i) {return z[i];};
  

  // optional member data accessors
  const double * gpstime_ptr() const {return gpstime;};
  const unsigned short * intensity_ptr() const {return intensity;};
  const unsigned char * classification_ptr() const {return classification;};
  const rgb48 * RGB_ptr() const {return RGB;};

  // user-defined member data accessors
  const double * data(std::string field) const {return extra_data.at(field);};

  //// old stuff
  void print_summary() const;
  void print_detailed() const;

  // mutators
  void add_intensity();
  void add_classification();
  void add_gpstime();
  void add_RGB();  
  void add_extra_data(std::string fieldname);

  PointCloud * subset(bool *keep);
  PointCloud * subset(unsigned int * keep_inds, unsigned int keep_count);
  PointCloud * copy();

  static PointCloud * read_LAS(std::string filename, unsigned int byte_offset=0);
  void write_LAS(std::string filename);


protected:

private:
  double *x, *y, *z, *gpstime;
  double xmin, xmax, ymin, ymax, zmin, zmax, gpst_min, gpst_max;
  unsigned int pointcount;
  unsigned short *intensity;
  unsigned char *classification;
  rgb48 * RGB;
  std::vector<std::string> extra_data_names;
  std::map<std::string, double *> extra_data;

  void calc_extents();  

  void read_LAS_internal(std::string filename, unsigned int byte_offset=0);

};

#pragma pack(push,1)
struct las_pt_0{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
};

struct las_pt_1{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
};

struct las_pt_2{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  unsigned short Red;
  unsigned short Green;
  unsigned short Blue;
};

struct las_pt_3{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
  unsigned short Red;
  unsigned short Green;
  unsigned short Blue;
};

struct las_pt_4{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
  unsigned char WaveP_Desc_Idx;
  unsigned long WaveP_Byte_Offset;
  unsigned int WaveP_Size;
  float WaveF_Location;
  float X_t;
  float Y_t;
  float Z_t;
};

struct las_pt_5{
  int X;
  int Y;
  int Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
  unsigned short Red;
  unsigned short Green;
  unsigned short Blue;
  unsigned char WaveP_Desc_Idx;
  unsigned long WaveP_Byte_Offset;
  unsigned int WaveP_Size;
  float WaveF_Location;
  float X_t;
  float Y_t;
  float Z_t;
};
#pragma pack(pop)

#endif
