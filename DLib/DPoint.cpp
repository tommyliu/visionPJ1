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
#include "DPoint.h"
#include <iostream>
#include <assert.h>
#include <stdlib.h>
using namespace std;

DPoint operator+(const DPoint &p1, const DPoint &p2)
{
  return DPoint(p1.row() + p2.row(), p1.col() + p2.col());
}

DPoint operator-(const DPoint &p1, const DPoint &p2)
{
  return DPoint(p1.row() - p2.row(), p1.col() - p2.col());
}

DPoint operator*(const DPoint &p1, const DPoint &p2)
{
  return DPoint(p1.row() * p2.row(), p1.col() * p2.col());
}

DPoint operator/(const DPoint &p1, int p2)
{
  return DPoint(p1.row() / p2, p1.col() / p2);
}


ostream &operator<<(ostream &os, const DPoint &p)
{
  os << "(" << p.row() << ", " << p.col() << ")";
  return os;
}

std::istream &operator>>(std::istream &is, DPoint &p)
{
  char c;

  do { is >> c; } while(c == ' ');

  assert(c=='(');

  is >> p._row;

  do { is >> c; } while(c == ' ');

  assert(c==',');

  is >> p._col;

  do { is >> c; } while(c == ' ');

  assert(c==')');

  return is;
}
