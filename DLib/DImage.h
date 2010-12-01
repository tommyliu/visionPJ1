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
#ifndef __IMAGE_H_
#define __IMAGE_H_

#include "DRect.h"
#include "DPlane.h"
#include "DIntegralImage.h"

class DImage 
{
  public:
  DImage() 
    {
      plane_count=0;
      plane_data=0;
    }
  
  DImage(int _planes, int _rows, int _cols)
    {
      plane_data = new DPlane *[_planes];
      plane_count = _planes;
      
      for(int i=0; i<_planes; i++)
	plane_data[i] = new DPlane(_rows, _cols);
    }

    DImage(int _planes)
    {
       plane_data = new DPlane *[_planes];
       plane_count = _planes;

       for(int i=0; i<_planes; i++)
          plane_data[i] = new DPlane;
    }

    ~DImage()
      {
	if(plane_data)
	  {
	    for(int i=0; i<planes(); i++)
	      delete plane_data[i];
	    
	    delete[] plane_data;
	  }
      }

    DImage(const DImage &other)
      {
	plane_data = 0;
	plane_count=0;

	*this = other;
      }

    DImage(const DPlane &other)
     {
       int _planes = 1;
       plane_data = new DPlane *[_planes];
       plane_count = _planes;

       plane_data[0] = new DPlane;
       (*plane_data[0]) = other;
     }

    DPlane &operator[](int plane) const
    {
       return *plane_data[plane];
    }

    DImage &operator=(const double D)
      {
	for(int i=0; i<plane_count; i++)
	  (*this)[i] = D;

	return *this;
      }

    DImage &operator=(const DImage &other)
      {
       if(plane_data)
	 {
	   for(int i=0; i<plane_count; i++)
	     delete plane_data[i];

	   delete[] plane_data;
	 }

	plane_count = other.planes();
	plane_data = new DPlane *[other.planes()];

	for(int i=0; i<plane_count; i++)
	    plane_data[i] = new DPlane(other[i]);

	return *this;

      }

    DPlane get_luma_plane() const
    {
       return (*this)[0]*0.2+(*this)[1]*0.7+(*this)[2]*0.1;
    }

    int rows() const { return plane_data[0]->rows(); }
    int cols() const { return plane_data[0]->cols(); }
    int planes() const { return plane_count; }

    friend std::istream &operator>>(std::istream &is, DImage &matrix);


  protected:
    int plane_count;
    DPlane **plane_data;

};


#endif
