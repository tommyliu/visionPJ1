//
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
#ifndef _DPOINT_H_
#define _DPOINT_H_

#include <iostream>

class DPoint
{
 public:

  DPoint(int __row, int __col) : _row(__row), _col(__col) {}
  DPoint() {}

  int row() const { return _row; }
  int col() const { return _col; }

  void row(int a) { _row = a; }
  void col(int a) { _col = a; }

  friend DPoint operator+(const DPoint &p1, const DPoint &p2);
  friend DPoint operator-(const DPoint &p1, const DPoint &p2);
  friend DPoint operator*(const DPoint &p1, const DPoint &p2);
  friend DPoint operator/(const DPoint &p1, const DPoint &p2);
  friend DPoint operator/(const DPoint &p1, int p2);

  friend std::ostream &operator<<(std::ostream &os, const DPoint &p);
  friend std::istream &operator>>(std::istream &is, DPoint &p);

 protected:
  int _row, _col;

};

#endif
