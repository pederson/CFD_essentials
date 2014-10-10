#ifndef _HULL_H
#define _HULL_H

#include <stdlib.h>
#include <vector>

class Point{
public:
  // constructor
  Point();

  // destructor
  ~Point();

  // data members
  double x, y, z;
  unsigned char ndims;
}

class Hull{
public:
  // constructor
  Hull();

  // destructor
  ~Hull();

  // data members
  vector<unsigned int> inds; // optional: indices from which the points came
  vector<Point> points; // ordered points that make up the hull

  // member functions
  construct_convex(vector<Point> allpts)
}
#endif
