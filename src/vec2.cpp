#include "vec2.hpp"

#include <algorithm> // for std::max
#include <stdlib.h>  // for abs


Vec2::Vec2() {}

Vec2::Vec2(int x, int y)
  : x_(x), y_(y) {}

Vec2 operator- (const Vec2& a, const Vec2& b) {
  return (Vec2){a.x_-b.x_, a.y_-b.y_};
}

bool operator== (const Vec2& a, const Vec2& b) {
  return (a.x_ == b.x_) && (a.y_ == b.y_);
}

bool operator!= (const Vec2& a, const Vec2& b) {
  return ! (a==b);
}

Vec2 Abs(const Vec2& a) {
  return (Vec2){abs(a.x_), abs(a.y_)};
}


double Max(const Vec2& a) {
  return std::max(a.x_, a.y_);
}

