#ifndef _HULL_H
#define _HULL_H

#include <stdlib.h>

#include <iostream>
#include <vector>
#include <stack>

class Point{
public:
  // constructor
  Point();

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

  // data members
  vector<unsigned int> inds; // optional: indices from which the points came
  vector<Point> points; // ordered points that make up the hull

  // member functions
  void construct_convex(vector<Point> allpts);
  void construct_convex(double * x, double * y, unsigned int numpoints);

protected:

private:
  Point p0; // used for sorting

  Point next_to_top(stack<Point> &S);
  int swap(Point &p1, Point &p2);
  int orientation(Point p, Point q, Point r);
  int compare(const void *vp1, const void *vp2);

};
#endif
