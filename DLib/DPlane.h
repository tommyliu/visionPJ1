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
#ifndef _DPLANE_H_
#define _DPLANE_H_

#include "DMatrix.h"
#include "DRect.h"
#include <string.h>

template<class T>
class _DPlane : public _DMatrix<T>
{

 public:
  _DPlane(const _DMatrix<T> &other)
    {
      _DMatrix<T>::_rows = other.rows();
      _DMatrix<T>::_cols = other.cols();
      
      _DMatrix<T>::data = 0;
      _DMatrix<T>::data_area = 0;
            
      _DMatrix<T>::initialize_storage();

      memcpy(_DMatrix<T>::data_area, other[0], 
	     _DMatrix<T>::_rows * _DMatrix<T>::_cols * sizeof(T));
    }
    
  _DPlane(int _rows, int _cols) : _DMatrix<T>(_rows, _cols)
    {
    }
  
  _DPlane() : _DMatrix<T>() {}
  
  void draw(const DRect &rect, T color);

  _DPlane<T> get_x_gradient() const;
  _DPlane<T> get_y_gradient() const;

  T find_min_val(DPoint &best_pt);
  T find_max_val(DPoint &best_pt);

  _DPlane<T> convolve(const _DPlane<T> &kernel) const;
  _DPlane<T> subsample(int row_factor, int col_factor);

  _DPlane<T> bilinear_interpolate(T row_offset, T col_offset) const;
  _DPlane<T> rotate_image(double angle, T bg_val=1e100) const;
  _DPlane<T> rotate_image_nn(double angle) const;

  _DPlane<T> rescale(const _DPlane<T> &img, T row_scale, T col_scale) const;
  
  _DPlane<T> binary_thin(void) const;

  _DMatrix<T> &operator=(const T &val);
};

typedef _DPlane<double> DPlane;
typedef _DPlane<short> DShortPlane;

                        



#endif
