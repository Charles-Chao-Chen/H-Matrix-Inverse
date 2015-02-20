#include "dim2.hpp"

#include <algorithm> // for std::max
#include <assert.h>
#include <stdlib.h>  // for abs


Dim2::Dim2() {}

Dim2::Dim2(int x, int y)
  : x_(x), y_(y) {}

int Dim2::area() const {
  return x_*y_;
}

Dim2 Dim2::x_bisection() const {
  int firstHalf  = x_ / 2;
  int secondHalf = x_ - firstHalf;
  return Dim2( firstHalf, secondHalf );
}

Dim2 Dim2::y_bisection() const {
  int firstHalf  = y_ / 2;
  int secondHalf = y_ - firstHalf;
  return Dim2( firstHalf, secondHalf );
}

int Dim2::operator[] (const int idx) {
  switch( idx ) {
  case 0:
    return x_;
    break;
  case 1:
    return y_;
    break;
  default:
    assert(false);
    break;
  }
}

Dim2 operator- (const Dim2& a, const Dim2& b) {
  return (Dim2){a.x_-b.x_, a.y_-b.y_};
}

Dim2 operator+ (const Dim2& a, const Dim2& b) {
  return (Dim2){a.x_+b.x_, a.y_+b.y_};
}

Dim2 operator* (const int scale, const Dim2& a) {
  return (Dim2){scale*a.x_, scale*a.y_};
}

Dim2 operator/ (const Dim2& a, const int scale) {
  return (Dim2){a.x_/scale, a.y_/scale};
}

bool operator== (const Dim2& a, const Dim2& b) {
  return (a.x_ == b.x_) && (a.y_ == b.y_);
}

bool operator!= (const Dim2& a, const Dim2& b) {
  return ! (a==b);
}

Dim2 Abs(const Dim2& a) {
  return (Dim2){abs(a.x_), abs(a.y_)};
}


double Max(const Dim2& a) {
  return std::max(a.x_, a.y_);
}

