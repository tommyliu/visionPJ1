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
#include "DRect.h"


std::ostream & operator<<(std::ostream &os, const DRect &rect)
{
  os << "[" << rect.top_left() << ", " << rect.bottom_right() << "]";
  return os;
}

std::istream & operator>>(std::istream &is, DRect &rect)
{
  char c;

  do { is >> c;  } while(c == ' ');

  assert(c=='[');

  is >> rect._top_left;

  do { is >> c; } while(c == ' ');

  assert(c==',');

  is >> rect._bottom_right;

  do { is >> c; } while(c == ' ');

  assert(c==']');
  
  return is;
}
/*
DRect compute_overlap(const DRect &r2)
{
    if(bottom() < r2.top() || r2.bottom() < top() ||
     right() < r2.left() || r2.right() < left())
    return 0;
  
  double row_overlap = min(r1.bottom() - r2.top(), r2.bottom() - r1.top());
  double col_overlap = min(r1.right() - r2.left(), r2.right() - r1.left());

  return row_overlap * col_overlap;
  

  return *this;
}*/
