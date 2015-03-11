#ifndef DIM2_H
#define DIM2_H

// grid points in 2d with integer coordinates
class Point2 {
public:
  Point2();
  Point2(int,int);

  // two static methods :
  //  (1) the point with absolution coordinates
  //  (2) the larger coordinates
  static Point2 abs(const Point2&);
  static int max(const Point2&);

public:
  int x_;
  int y_;
};

/* --- operator overloading --- */

Point2 operator- (const Point2& a, const Point2& b);

Point2 operator+ (const Point2& a, const Point2& b);

Point2 operator* (const int scale, const Point2& a);

Point2 operator/ ( const Point2& a, const int scale);

bool operator== (const Point2& a, const Point2& b);

bool operator!= (const Point2& a, const Point2& b);

// rectangle in 2 dimensions with width and height in integers
class Rect2 {
public:
  // constructor
  Rect2();
  Rect2(int,int);

  // member functions

  // return the area of the rectangle
  int area() const;

  // the first and second half of the width
  Rect2 x_bisection() const;

  // the first and second half of the height
  Rect2 y_bisection() const;
  
  /* --- operator overloading --- */
  int operator[] (const int);
  
private:
  int x[2];
};

#endif // DIM2_H
