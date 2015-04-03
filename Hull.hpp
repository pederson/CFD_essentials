#ifndef _HULL_H
#define _HULL_H

#include <stdlib.h>

#include <iostream>
#include <vector>
#include <stack>

enum HullType{CONVEX=0, CONCAVE=1};
enum Orient{CW=0, CCW=1};

class Point{
  friend class Hull;
public:
  // constructor
  Point(double _x=0.0, double _y=0.0, double _z=0.0)
    :x(_x),
     y(_y),
     z(_z)
     {_ndims = 3;}

  // destructor
  ~Point();

  // inspectors
  //double x() const {return x;};
  //double y() const {return y;};
  //double z() const {return z;};
  void print_summary() const;

  // mutators


  static int dist(const Point & p1, const Point & p2);

  

private:
  // data members
  double x, y, z;
  unsigned char _ndims;

};

class Hull{
public:
  // constructor
  Hull();
  Hull(std::vector<Point> ordered_points);
  Hull(const double * x, const double * y, unsigned int numpoints);

  // destructor
  ~Hull();

  // inspectors
  void print_summary() const;
  void print_detailed() const;
  bool contains_point(Point query) const;

  // mutators
  void construct_convex(std::vector<Point> allpts);
  void construct_convex(const double * x, const double * y, unsigned int numpoints);

  

  Hull * Union(Hull * other_Hull);
  Hull * Intersect(Hull * other_Hull);
  Hull * Stamper(std::vector<Point>); 

private:

  // data members
  std::vector<unsigned int> inds; // optional: indices from which the points came
  std::vector<Point> hull_points; // ordered points that make up the hull
  HullType type;
  Orient direction; // direction of ordering (CCW or CW)
  double xmin, xmax, ymin, ymax;

  void calc_extents();

  static Point next_to_top(std::stack<Point> &S);
  static int swap(Point &p1, Point &p2);
  static int orientation(Point p, Point q, Point r);
  static int compare(const void *vp1, const void *vp2);
  static bool lines_intersect_query(Point p1, Point q1, Point p2, Point q2);
  static bool on_segment_query(Point p, Point q, Point r);

};

extern Point p0; // global used for sorting

#endif
