#include "vec2.hpp"

#include <algorithm> // for std::max
#include <stdlib.h>  // for abs


Vec2 operator- (const Vec2& a, const Vec2& b) {
  return (Vec2){a.x-b.x, a.y-b.y};
}

bool operator== (const Vec2& a, const Vec2& b) {
  return (a.x == b.x) && (a.y == b.y);
}

bool operator!= (const Vec2& a, const Vec2& b) {
  return ! (a==b);
}

Vec2 Abs(const Vec2& a) {
  return (Vec2){abs(a.x), abs(a.y)};
}


double Max(const Vec2& a) {
  return std::max(a.x, a.y);
}

