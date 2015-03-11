#include "dim2.hpp"

#include <algorithm> // for std::max and std::abs
#include <assert.h>

Point2::Point2()
  : x_(0), y_(0) {}

Point2::Point2(int x, int y)
  : x_(x), y_(y) {}

/*static*/ Point2 Point2::abs(const Point2& a) {
  return Point2( std::abs(a.x_), std::abs(a.y_) );
}

/*static*/ int Point2::max(const Point2& a) {
  return std::max(a.x_, a.y_);
}

Point2 operator- (const Point2& a, const Point2& b) {
  return Point2(a.x_-b.x_, a.y_-b.y_);
}

Point2 operator+ (const Point2& a, const Point2& b) {
  return Point2(a.x_+b.x_, a.y_+b.y_);
}

Point2 operator* (const int scale, const Point2& a) {
  return Point2(scale*a.x_, scale*a.y_);
}

Point2 operator/ (const Point2& a, const int scale) {
  return Point2(a.x_/scale, a.y_/scale);
}

bool operator== (const Point2& a, const Point2& b) {
  return (a.x_ == b.x_) && (a.y_ == b.y_);
}

bool operator!= (const Point2& a, const Point2& b) {
  return ! (a==b);
}

Rect2::Rect2() { x[0]=x[1]=0; }

Rect2::Rect2(int a, int b) { x[0]=a, x[1]=b;}

int Rect2::area() const {
  return x[0]*x[1];
}

Rect2 Rect2::x_bisection() const {
  int firstHalf  = x[0] / 2;
  int secondHalf = x[0] - firstHalf;
  return Rect2( firstHalf, secondHalf );
}

Rect2 Rect2::y_bisection() const {
  int firstHalf  = x[1] / 2;
  int secondHalf = x[1] - firstHalf;
  return Rect2( firstHalf, secondHalf );
}

int Rect2::operator[] (const int idx) {
#ifdef DEBUG
  assert( idx==0 || idx==1 );
#endif
  return x[idx];
}

