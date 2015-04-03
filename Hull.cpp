#include "Hull.hpp"

//#define _TEST_

using namespace std;

Point p0;


// destructor
Point::~Point(){

}

int Point::dist(const Point & p1, const Point & p2){
  return (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z);
}

void Point::print_summary() const{
  cout << "x: " << x << ", y: " << y << ", z: " << z << endl;
  return;
}

Hull::Hull(){

}

Hull::Hull(std::vector<Point> ordered_points){

  if (ordered_points.size() < 3){
    cout << "ERROR: need more than 3 points for a hull" << endl;
    throw -1;
  }

  type = CONCAVE;

  hull_points = ordered_points;

  // check the direction
  direction = CCW;

  // calculate the max and min
  calc_extents();
}

Hull::Hull(const double * x, const double * y, unsigned int numpoints){
  cout << "OH NOES THIS CTOR ISNT YET IMPLEMENTED" << endl;
}

  // destructor
Hull::~Hull(){

}

void Hull::print_summary() const{
  cout << " " << endl;
  cout << "************** Hull summary **************" << endl;
  if (hull_points.size() < 1) {
    cout << "Hull is empty!" << endl;
    return;
  }

  cout << "  type: ";
  if (type==CONVEX) cout << "CONVEX" << endl;
  else if (type==CONCAVE) cout << "CONCAVE" << endl;
  else cout << "UNKNOWN" << endl;
  cout << "  pointcount: " << hull_points.size() << endl;
  cout << "  x extents: [" << xmin << ", " << xmax << "]" << endl;
  cout << "  y extents: [" << ymin << ", " << ymax << "]" << endl;
  cout << "******************************************" << endl;
  cout << " " << endl;

  return;
}

void Hull::print_detailed() const{
  cout << " " << endl;
  cout << "************** Hull details **************" << endl;
  if (hull_points.size() < 1) {
    cout << "Hull is empty!" << endl;
    return;
  }

  cout << "  type: ";
  if (type==CONVEX) cout << "CONVEX" << endl;
  else if (type==CONCAVE) cout << "CONCAVE" << endl;
  else cout << "UNKNOWN" << endl;
  cout << "  pointcount: " << hull_points.size() << endl;
  cout << "  x extents: [" << xmin << ", " << xmax << "]" << endl;
  cout << "  y extents: [" << ymin << ", " << ymax << "]" << endl;
  cout << "******************************************" << endl;
  cout << " " << endl;

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
  //p0.print_summary();
  //for (int i=0; i<allpts.size(); i++) allpts[i].print_summary();
  //cout << "about to enter qsort" << sizeof(Point) << endl;
  qsort(&allpts[1], allpts.size()-1, sizeof(Point), Hull::compare);

  //cout << "finished qsort" << endl;
  //for (int i=0; i<allpts.size(); i++) allpts[i].print_summary();
  //create an empty stack and push first three points to it
  stack<Point> S;
  S.push(allpts[0]);
  S.push(allpts[1]);
  S.push(allpts[2]);

  //cout << "about to process remaining" << endl;
  int orient_val;
  // process remaining n-3 points
  for (int i=3; i<allsize; i++){
    // keep removing top while the angle formed by points next to top, 
    // top, and points[i] makes a non-left turn
    //cout << "i: " << i << endl;

    orient_val = orientation(next_to_top(S), S.top(), allpts[i]);
    while (orient_val != 2 ){
      if (orient_val == 0) break;
      //cout << "size of S: " << S.size() << endl;
      //cout << "gonna pop from S" << endl;
      S.pop();
      //S.top().print_summary();
      orient_val = orientation(next_to_top(S), S.top(), allpts[i]);
    } 
    //if (orient_val == 0) S.pop();

    //cout << "made it to the push" << endl;
    S.push(allpts[i]);
  }

  //cout << "finished processing stack has " << S.size() << " points" << endl;
  
  // now stack has the output points...move it into the member data
  hull_points.resize(S.size());
  int numstack = S.size();
  for (int i=0; i<numstack; i++){
    hull_points[i] = S.top();
    S.pop();
  }
  type = CONVEX;
  direction = CW;

  calc_extents();

  return;
}

void Hull::construct_convex(const double * x, const double * y, unsigned int numpoints){
  // declare vars
  vector<Point> setpts;

  setpts.resize(numpoints);
  // create a vector of points from the x,y pairs
  for (unsigned int i=0; i<numpoints; i++){
    setpts[i] = Point(x[i], y[i]);
  }

  // send it to the core function
  construct_convex(setpts);

  return;
}

bool Hull::contains_point(Point p) const{
    // There must be at least 3 vertices in polygon[]
    if (hull_points.size() < 3){
      cout << "WARNING: hull is less than 3 points" << endl;
      return false;
    }

    // Check if the point is outside the min/max
    if (p.x > xmax || p.x < xmin || p.y > ymax || p.y < ymin){
      return false;
    }
 
    // Create a point for line segment from p to infinite
    Point extreme = {xmax+100000, p.y};
 
    // Count intersections of the above line with sides of hull
    unsigned int count = 0, i = 0;
    do
    {
        unsigned int next = (i+1)%hull_points.size();
 
        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'hull[i]' to 'hull[next]'
        if (lines_intersect_query(hull_points.at(i), hull_points.at(next), p, extreme))
        {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(hull_points.at(i), p, hull_points.at(next)) == 0)
               return on_segment_query(hull_points.at(i), p, hull_points.at(next));
 
            count++;
        }
        i = next;
    } while (i != 0);
 
    // Return true if count is odd, false otherwise
    return (count%2 ==1);  // Same as (count%2 == 1) (count&1)
}

Point Hull::next_to_top(stack<Point> &S){
  //cout << "about to do next to top" ;
  Point p = S.top();
  S.pop();
  Point res = S.top();
  S.push(p);
  //cout << "next to top succeeded" << endl;
  return res;
}

int Hull::swap(Point &p1, Point &p2){
  Point temp = p1;
  p1 = p2;
  p2 = temp;
}

int Hull::orientation(Point p, Point q, Point r){
  double val = (q.y - p.y)*(r.x - q.x) - (q.x - p.x)*(r.y - q.y);
  
  /*
  p.print_summary();
  q.print_summary();
  r.print_summary();
  cout << "orientation succeeded... val is " << val << endl;
  */
  

  if (val == 0) return 0; // collinear
  return (val > 0)? 1: 2; // 1=clockwise (right turn) or 2=counterclockwise
}

int Hull::compare(const void *vp1, const void *vp2){
  Point *p1 = (Point *) vp1;
  Point *p2 = (Point *) vp2;

  // find orientation
  int o = orientation(p0, *p1, *p2);
  if (p0.x == p1->x && p0.x == p2->x) return (Point::dist(p0, *p2) >= Point::dist(p0, *p1))? 1 : -1;
  if (o == 0) return (Point::dist(p0, *p2) >= Point::dist(p0, *p1))? -1 : 1;

  return (o == 2)? -1: 1;
}

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool Hull::on_segment_query(Point p, Point q, Point r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
            q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
        return true;
    return false;
}

bool Hull::lines_intersect_query(Point p1, Point q1, Point p2, Point q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
 
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && on_segment_query(p1, p2, q1)) return true;
 
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && on_segment_query(p1, q2, q1)) return true;
 
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && on_segment_query(p2, p1, q2)) return true;
 
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && on_segment_query(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}


#ifdef _TEST_


int main(int argc, char * argv[]){
  // declare vars
  double *setx, *sety;
  unsigned int N=31;
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

  cout << "about to construct hull" << endl;

  // create a convex hull
  testhull->construct_convex(setx, sety, N*N);

  cout << "constructed the hull" << endl;

  // output summary
  testhull->print_summary();

  // test the intersect function
  cout << "hull contains (0.5, 0.5)? " << (testhull->contains_point({0.5, 0.5})? "yup!" : "no :X") << endl;
  cout << "hull contains (-1.0, 0.5)? " << (testhull->contains_point({-1.0, 0.5})? "yup!" : "no :X") << endl;
  cout << "hull contains (0, 0)? " << (testhull->contains_point({0.0, 0.0})? "yup!" : "no :X") << endl;

  //for (unsigned int i=0; i<testhull->hull_points.size(); i++) testhull->hull_points[i].print_summary();

  delete testhull;

  return 0;
}
#endif