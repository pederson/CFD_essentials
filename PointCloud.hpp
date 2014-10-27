#ifndef _POINTCLOUD_H
#define _POINTCLOUD_H

#include <fstream>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>

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

  PointCloud(unsigned int numpts=0);
  ~PointCloud();

  double *x, *y, *z, *gpstime;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  unsigned int pointcount;
  unsigned short *intensity;
  unsigned char *classification;
  rgb48 * RGB;

  void print_summary();
  void calc_extents();

  void add_intensity();
  void add_classification();
  void add_gpstime();
  void add_RGB();  

  PointCloud * subset(bool *keep);
  PointCloud * subset(unsigned int * keep_inds, unsigned int keep_count);

  static PointCloud * read_LAS(char * filename, unsigned int byte_offset=0);

protected:

private:

  void write_LAS(char * filename);

};

#pragma pack(push,1)
struct las_pt_0{
  long X;
  long Y;
  long Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
};

struct las_pt_1{
  long X;
  long Y;
  long Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
};

struct las_pt_2{
  long X;
  long Y;
  long Z;
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
  long X;
  long Y;
  long Z;
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
  long X;
  long Y;
  long Z;
  unsigned short Intensity;
  unsigned char Return_Info;
  unsigned char Classification;
  char Scan_Angle_Rank;
  unsigned char User_Data;
  unsigned short Point_Source_ID;
  double GPSTime;
  unsigned char WaveP_Desc_Idx;
  unsigned long long WaveP_Byte_Offset;
  unsigned long WaveP_Size;
  float WaveF_Location;
  float X_t;
  float Y_t;
  float Z_t;
};

struct las_pt_5{
  long X;
  long Y;
  long Z;
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
  unsigned long long WaveP_Byte_Offset;
  unsigned long WaveP_Size;
  float WaveF_Location;
  float X_t;
  float Y_t;
  float Z_t;
};
#pragma pack(pop)

#endif
