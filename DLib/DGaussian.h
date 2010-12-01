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
#include <DPlane.h>

class DGaussian : public DPlane
{
  public:
  DGaussian(int width, double sigma) : DPlane(1, width)
  {
    assert(width % 2 == 1);
    double norm = 1.0/(sigma * sqrt(2*M_PI));
    int half_width = width/2;
    for(int i=0; i<width; i++)
      (*this)[0][i] = norm * exp(-pow((i-half_width)-sigma,2)/(2*sigma*sigma));
  }

  DGaussian(int height, int width, double sigma) : DPlane(height, width)
  {
    assert(width % 2 == 1 && height % 2 == 1);
    double norm = 1.0/(sigma * sqrt(2*M_PI));
    int half_width = width/2;
    int half_height = height/2;
    for(int i=0; i<height; i++)
      for(int j=0; j<width; j++)
        (*this)[i][j] = norm * norm *
                      exp(-(pow(i-half_height,2) + pow(j-half_width, 2))/
                      (2*sigma*sigma));
  }


}; 
