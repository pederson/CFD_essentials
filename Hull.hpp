#ifndef _HULL_H
#define _HULL_H

#include <stdlib.h>

#include <iostream>
#include <vector>
#include <stack>

enum HullType{CONVEX=0, CONCAVE=1};

class Point{
  friend class Hull;
public:
  // constructor
  Point();
  Point(double _x, double _y);
  Point(double _x, double _y, double _z);

  // destructor
  ~Point();

  static int dist(Point p1, Point p2);

private:
  // data members
  double x, y, z;
  unsigned char ndims;

};

class Hull{
public:
  // constructor
  Hull();

  // destructor
  ~Hull();

  void print_summary();
  void calc_extents();

  // data members
  std::vector<unsigned int> inds; // optional: indices from which the points came
  std::vector<Point> hull_points; // ordered points that make up the hull
  HullType type;
  double xmin, xmax, ymin, ymax;

  // member functions
  void construct_convex(std::vector<Point> allpts);
  void construct_convex(double * x, double * y, unsigned int numpoints);

  void printSummary();



protected:

private:

  static Point next_to_top(std::stack<Point> &S);
  static int swap(Point &p1, Point &p2);
  static int orientation(Point p, Point q, Point r);
  static int compare(const void *vp1, const void *vp2);

};

Point p0; // global used for sorting

#endif
