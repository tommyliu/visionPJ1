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
#ifndef _POINT_H_
#define _POINT_H_

#include "DPoint.h"
#include <iostream>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

class DRect
{
 public:
  DRect(const DPoint &__top_left, const DPoint &__bottom_right) : 
    _top_left(__top_left), _bottom_right(__bottom_right) {}
  DRect(int __top, int __left, int __bottom, int __right) :
    _top_left(DPoint(__top, __left)), _bottom_right(DPoint(__bottom, __right)) {}
  DRect()
    {
      _top_left=DPoint(-1,-1);
      _bottom_right=DPoint(-1,-1);
    }

  DPoint top_left() const { return _top_left; }
  DPoint bottom_right() const { return _bottom_right; }
  DPoint top_right() const { return DPoint(top(), right()); }
  DPoint bottom_left() const { return DPoint(bottom(), left()); }

  DPoint center_point() const 
  {
    return top_left() + size()/2;
  }

  int top() const { return _top_left.row(); }
  int bottom() const { return _bottom_right.row(); }
  int left() const { return _top_left.col(); }
  int right() const { return _bottom_right.col(); }

  void top(int a) { _top_left.row(a); }
  void bottom(int a) { _bottom_right.row(a); }
  void left(int a) { _top_left.col(a); }
  void right(int a) { _bottom_right.col(a); }

  int height() const { return bottom() - top() + 1; }
  int width() const { return right() - left() + 1; }
  DPoint size() const { return DPoint(height(), width()); }

  int area() const { return height() * width(); }

  //  DRect compute_overlap(const DRect &r2);

  DRect operator+(const DPoint &pt)
  {
    return DRect(_top_left + pt, _bottom_right + pt);
  }

  friend std::ostream & operator<<(std::ostream &os, const DRect &rect);
  friend std::istream & operator>>(std::istream &is, DRect &rect);

 protected:
  DPoint _top_left, _bottom_right;
};


#endif 
