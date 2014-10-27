#include "Hull.hpp"

#define _TEST_

using namespace std;

Point::Point(){

}

Point::Point(double _x, double _y){
  x = _x;
  y = _y;
  z = 0.0;
  ndims = 2;
}

Point::Point(double _x, double _y, double _z){
  x = _x;
  y = _y;
  z = _z;
  ndims = 3;
}

// destructor
Point::~Point(){

}

int Point::dist(Point p1, Point p2){
  return (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z);
}

Hull::Hull(){

}

  // destructor
Hull::~Hull(){

}

void Hull::print_summary(){
  if (hull_points.size() < 1) {
    cout << "Hull is empty!" << endl;
    return;
  }

  cout << "Hull summary:" << endl;
  cout << "  type: ";
  if (type==CONVEX) cout << "CONVEX" << endl;
  else if (type==CONCAVE) cout << "CONCAVE" << endl;
  else cout << "UNKNOWN" << endl;
  cout << "  pointcount: " << hull_points.size() << endl;
  cout << "  x extents: [" << xmin << ", " << xmax << "]" << endl;
  cout << "  y extents: [" << ymin << ", " << ymax << "]" << endl;

  return;
}

void Hull::calc_extents(){
  if (hull_points.size() < 1) {
    cout << "Hull is empty!" << endl;
    return;
  }

  xmin = hull_points[0].x; xmax = hull_points[0].x; ymin = hull_points[0].y; ymax = hull_points[0].y;
  for (unsigned int i=1; i<hull_points.size(); i++){
    if (hull_points[i].x < xmin) xmin = hull_points[i].x;
    if (hull_points[i].x > xmax) xmax = hull_points[i].x;
    if (hull_points[i].y < ymin) ymin = hull_points[i].y;
    if (hull_points[i].y > ymax) ymax = hull_points[i].y;
  }

  return;
}

// member functions
void Hull::construct_convex(vector<Point> allpts){

  unsigned int allsize = allpts.size();

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
  cout << "about to enter qsort" << endl;
  qsort(&allpts[1], allpts.size()-1, sizeof(Point), Hull::compare);

  cout << "finished qsort" << endl;
  //create an empty stack and push first three points to it
  stack<Point> S;
  S.push(allpts[0]);
  S.push(allpts[1]);
  S.push(allpts[2]);

  cout << "about to process remaining" << endl;

  // process remaining n-3 points
  for (int i=3; i<allsize; i++){
    // keep removing top while the angle formed by points next to top, 
    // top, and points[i] makes a non-left turn
    cout << "i: " << i << endl;

    while (orientation(next_to_top(S), S.top(), allpts[i]) != 2) S.pop();

    cout << "made it to the push" << endl;
    S.push(allpts[i]);
  }
  cout << "finished processing" << endl;

  // now stack has the output points...move it into the member data
  hull_points.resize(S.size());
  for (int i=0; i<S.size(); i++){
    hull_points[i] = S.top();
    S.pop();
  }
  type = CONVEX;

  calc_extents();

  return;
}

void Hull::construct_convex(double * x, double * y, unsigned int numpoints){
  // declare vars
  vector<Point> setpts;

  setpts.resize(numpoints);
  // create a vector of points from the x,y pairs
  for (unsigned int i=0; i<numpoints; i++){
    setpts[i] = Point(x[i], y[i]);
  }

  // send it to the core function

  cout << "about to enter core function " << endl;

  construct_convex(setpts);

  return;
}

Point Hull::next_to_top(stack<Point> &S){
  Point p = S.top();
  S.pop();
  Point res = S.top();
  S.push(p);
  cout << "next to top succeeded" << endl;
  return res;
}

int Hull::swap(Point &p1, Point &p2){
  Point temp = p1;
  p1 = p2;
  p2 = temp;
}

int Hull::orientation(Point p, Point q, Point r){
  int val = (q.y - p.y)*(r.x - q.x) - (q.x - p.x)*(r.y - q.y);

  cout << "orientation succeeded" << endl;

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


#ifdef _TEST_


int main(int argc, char * argv[]){
  // declare vars
  double *setx, *sety;
  unsigned int N=20;
  double space=1.0/double(N);
  Hull * testhull = new Hull();

  setx = new double[N*N];
  sety = new double[N*N];

  // create a bunch of values from (0,0) to (1,1)
  for (unsigned int i=0; i<N; i++){
    for (unsigned int j=0; j<N; j++){
      setx[i*N + j] = i*space;
      sety[i*N + j] = j*space;
    }
  }

  cout << "got here" << endl;

  // create a convex hull
  testhull->construct_convex(setx, sety, N*N);

  cout << "constructed the hull" << endl;


  // output summary
  testhull->print_summary();

  delete testhull;

  return 0;
}
#endif