#include "Hull.hpp"

using namespace std;

Point::Point(){

}

// destructor
Point::~Point(){

}

int Point::dist(Point p1, Point p2){
  return (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z);
}

Hull::Hull();

  // destructor
Hull::~Hull();

// member functions
void Hull::construct_convex(vector<Point> allpts){
  // find the bottommost point
  double ymin = allpts[0].y, min = 0;
  for (int i=1; i<allpts.size(); i++){
    double y = allpts[i].y;

    // pick the bottom-most or choose the left most point in case of tie
    if ((y < ymin) || (ymin == y && allpts[i].x < allpts[min].x)){
      ymin = allpts[i].y;
      min = i;
    }
  }

  // place the bottom-most point at first position
  swap(allpts[0], allpts[min]);

  // sort n-1 points with respect to the first point. A point p1 comes
  // before p2 in sorted output if p2 has larger polar angle (in CCW
  // direction) than p1
  p0 = allpts[0];
  qsort(&allpts[1], allpts.size()-1, sizeof(Point), compare);

  //create an empty stack and push first three points to it
  stack<Point> S;
  S.push(allpts[0]);
  S.push(allpts[1]);
  S.push(allpts[2]);

  // process remaining n-3 points
  for (int i=3; i<allpts.size(); i++){
    // keep removing top while the angle formed by points next to top, 
    // top, and points[i] makes a non-left turn
    while (orientation(next_to_top(S), S.top(), allpts[i]) != 2) S.pop();
    S.push(allpts[i]);
  }

  // now stack has the output points...move it into the member data
  hull_points.resize(S.size());
  for (int i=0; i<S.size(); i++){
    hull_points[i] = S.top();
    S.pop();
  }
}

void Hull::construct_convex(double * x, double * y, unsigned int numpoints);

Point Hull::next_to_top(stack<Point> &S){
  Point p = S.top();
  S.pop();
  Point res = S.top();
  S.push(p);
  return res;
}

int Hull::swap(Point &p1, Point &p2){
  Point temp = p1;
  p1 = p2;
  p2 = temp;
}

int Hull::orientation(Point p, Point q, Point r){
  int val = (q.y - p.y)*(r.x - q.x) - (q.x - p.x)*(r.y - q.y);

  if (val == 0) return 0; // collinear
  return (val > 0)? 1: 2; // clockwise or counterclockwise
}

int Hull::compare(const void *vp1, const void *vp2){
  Point *p1 = (Point *) vp1;
  Point *p2 = (Point *) vp2;

  // find orientation
  int o = orientation(p0, *p1, *p2);
  if (o == 0) return (Point::dist(p0, *p2) >= Point::dist(p0, *p1))? -1 : 1;

  return (o == 2)? -1: 1;
}


