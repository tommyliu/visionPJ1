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
#include <DImage.h>

std::istream &operator>>(std::istream &is, DImage &matrix)
{
  int planes;
  is >> planes;

  std::cout << "reading " << planes << std::endl;
  matrix = DImage(planes);

  for(int i=0; i<planes; i++)
    is >> matrix[i];

  return is;
}
