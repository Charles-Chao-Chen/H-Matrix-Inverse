#ifndef VEC2_H
#define VEC2_H

// Vec2 structure for variables in 2d
struct Vec2 {
  Vec2();
  Vec2(int,int);
  int x_;
  int y_;
};

/* --- operator overloading --- */

Vec2 operator- (const Vec2& a, const Vec2& b);

bool operator== (const Vec2& a, const Vec2& b);

bool operator!= (const Vec2& a, const Vec2& b);


Vec2 Abs(const Vec2&);

// return the max coordinate
double Max(const Vec2&);



#endif // VEC2_H
