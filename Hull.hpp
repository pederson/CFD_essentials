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
  Point();
  Point(double _x, double _y);
  Point(double _x, double _y, double _z);

  // destructor
  ~Point();

  static int dist(Point p1, Point p2);

  void print();

private:
  // data members
  double x, y, z;
  unsigned char ndims;

};

class Hull{
public:
  // constructor
  Hull();
  Hull(std::vector<Point> ordered_points);
  Hull(double * x, double * y, unsigned int numpoints);

  // destructor
  ~Hull();

  // member functions
  void print_summary() const;
  void print_detailed() const;

  void construct_convex(std::vector<Point> allpts);
  void construct_convex(double * x, double * y, unsigned int numpoints);

  bool contains_point(Point query) const;

  Hull * Union(Hull * other_Hull);
  Hull * Intersect(Hull * other_Hull);
  Hull * Stamper(std::vector<Point>); 



protected:

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
