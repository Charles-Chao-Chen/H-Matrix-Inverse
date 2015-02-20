#ifndef DIM2_H
#define DIM2_H

// for grid pints in 2d
//  note the coordinates are both integers
struct Dim2 {

  // constructor
  Dim2();
  Dim2(int,int);

  // member functions

  int area() const;

  Dim2 x_bisection() const;

  Dim2 y_bisection() const;
  
  /* --- operator overloading --- */

  // Dim2 object can be used as int[2]
  //  where [0] and [1] correspond to x_ and y_
  int  operator[] (const int);

  // member variables
  int x_;
  int y_;
};

/* --- operator overloading --- */

Dim2 operator- (const Dim2& a, const Dim2& b);

Dim2 operator+ (const Dim2& a, const Dim2& b);

Dim2 operator* (const int scale, const Dim2& a);

Dim2 operator/ ( const Dim2& a, const int scale);

bool operator== (const Dim2& a, const Dim2& b);

bool operator!= (const Dim2& a, const Dim2& b);


Dim2 Abs(const Dim2&);

// return the max coordinate
double Max(const Dim2&);



#endif // DIM2_H
