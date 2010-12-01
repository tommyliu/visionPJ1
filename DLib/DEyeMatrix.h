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
#include <DMatrix.h>

class DEyeMatrix : public DMatrix
{
  public:
    DEyeMatrix(int _rows, int _cols) : DMatrix(_rows, _cols)
    {
      for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++)
        {
          if(i == j)
            (*this)[i][j] = 1;
          else
            (*this)[i][j] = 0;
        }
    }


};
