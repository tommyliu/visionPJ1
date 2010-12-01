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
#include "DPlane.h"
#include "thin_lookup.h"

template<class T>
void _DPlane<T>::draw(const DRect &rect, T color)
{
  for(int i=rect.top(); i <= rect.bottom(); i++)
    {
      if(i >= 0 && i < _DMatrix<T>::rows() && rect.left() >= 0 && rect.left() < _DMatrix<T>::cols())
	 (*this)[i][rect.left()] = color;
      if(i >= 0 && i < _DMatrix<T>::rows() && rect.right() >= 0 && rect.right() < _DMatrix<T>::cols())
	(*this)[i][rect.right()] = color;
    }

  for(int j=rect.left(); j <= rect.right(); j++)
    {
      if(rect.top() >= 0 && rect.top() < _DMatrix<T>::rows() && j >= 0 && j < _DMatrix<T>::cols())
	(*this)[rect.top()][j] = color;
      if(rect.bottom() >= 0 && rect.bottom() < _DMatrix<T>::rows() && j >= 0 && j < _DMatrix<T>::cols())
	(*this)[rect.bottom()][j] = color;
    }
}

template<class T>
_DMatrix<T> &_DPlane<T>::operator=(const T &val)
{
  T *cp = _DMatrix<T>::data_area;
  for(int i=0; i<_DMatrix<T>::rows(); i++)
    for(int j=0; j<_DMatrix<T>::cols(); j++, cp++)
      *cp = val;

  return (*this);
}

template<class T>
T _DPlane<T>::find_max_val(DPoint &best_pt)
{
  int _rows = _DMatrix<T>::rows(), _cols = _DMatrix<T>::cols();
  T *cp = (*this)[0];
  T max_val = *cp;
  for(int i=0; i < _rows; i++)
    for(int j=0; j < _cols; j++, cp++)
      if(*cp > max_val)
	{
	  max_val = *cp;
	  best_pt = DPoint(i, j);
	}
    
  return max_val;
}

template<class T>
T _DPlane<T>::find_min_val(DPoint &best_pt)
{
  int _rows = _DMatrix<T>::rows(), _cols = _DMatrix<T>::cols();
  T *cp = (*this)[0];
  T min_val = *cp;
  for(int i=0; i < _rows; i++)
    for(int j=0; j < _cols; j++, cp++)
      if(*cp < min_val)
	{
	  min_val = *cp;
	  best_pt = DPoint(i, j);
	}
    
  return min_val;
}

template<class T>
_DPlane<T> _DPlane<T>::get_x_gradient() const
{
  _DPlane<T> result(_DMatrix<T>::rows()-1, _DMatrix<T>::cols()-1);

  result=0;

  for(int i=0; i<_DMatrix<T>::rows()-1; i++)
    for(int j=0; j<_DMatrix<T>::cols()-1; j++)
      result[i][j] = (*this)[i+1][j+1] - (*this)[i][j];

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::get_y_gradient() const
{
  _DPlane<T> result(_DMatrix<T>::rows()-1, _DMatrix<T>::cols()-1);

  result=0;

  for(int i=0; i<_DMatrix<T>::rows()-1; i++)
    for(int j=0; j<_DMatrix<T>::cols()-1; j++)
      result[i][j] = (*this)[i+1][j] - (*this)[i][j+1];

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::convolve(const _DPlane<T> &kernel) const
{
  _DPlane<T> result(_DMatrix<T>::rows(), _DMatrix<T>::cols());

  int half_height = kernel._DMatrix<T>::rows()/2;
  int half_width = kernel._DMatrix<T>::cols()/2;

  for(int i=0; i<_DMatrix<T>::rows(); i++)
    for(int j=0; j<_DMatrix<T>::cols(); j++)
    {
      T sum=0, total=0;
      for(int k=-half_height; k<=half_height; k++)
        for(int l=-half_width; l<=half_width; l++)
          if(i+k >=0 && j+l >=0 && 
             i+k < _DMatrix<T>::rows() && j+l < _DMatrix<T>::cols()) 
          {
            sum += kernel[k+half_height][l+half_width] * (*this)[i+k][j+l];
            total += kernel[k+half_height][l+half_width];
          }
      result[i][j] = sum/total;
    }

  return result;
}

template<class T>
_DPlane<T> _DPlane<T>::subsample(int row_factor, int col_factor)
{
  int r = (int) ceil(T(_DMatrix<T>::rows()) / row_factor);
  int c = (int) ceil(T(_DMatrix<T>::cols()) / col_factor);
  _DPlane<T> result(r, c);

  T *out_cp = result[0];
  int n=0;
  for(int i=0; i<_DMatrix<T>::rows(); i+=row_factor)
    for(int j=0; j<_DMatrix<T>::cols(); j+=col_factor, out_cp++,n++)
      *out_cp = (*this)[i][j];

  return result;
}

double float_part(double i)
{
  return i - floor(i);
}

template<class T>
_DPlane<T> _DPlane<T>::bilinear_interpolate(T row_offset, T col_offset) const
{
  assert(row_offset >= 0 && col_offset >= 0);
  assert(row_offset < 1 && col_offset < 1);

  int new_rows = _DMatrix<T>::rows() - (int)ceil(row_offset);
  int new_cols = _DMatrix<T>::cols() - (int)ceil(col_offset);

  _DPlane<T> result(new_rows, new_cols);

  T row_pos, col_pos;

  row_pos = row_offset;
  for(int i=0; i<new_rows; i++, row_pos++)
    {
      col_pos = col_offset;
      for(int j=0; j<new_cols; j++, col_pos++)
	{
	  double row_fraction = float_part(row_pos);
	  double col_fraction = float_part(col_pos);
	  int row_int = (int)floor(row_pos);
	  int col_int = (int)floor(col_pos);
	  
	  result[i][j] = 
	    T((1.0-row_fraction) * (1.0-col_fraction) * (*this)[row_int][col_int] +
	    (1.0-row_fraction) * (    col_fraction) * (*this)[row_int][col_int+1] +
	    (    row_fraction) * (1.0-col_fraction) * (*this)[row_int+1][col_int] +
	      (    row_fraction) * (    col_fraction) * (*this)[row_int+1][col_int+1]);
	    
	}
    }

  return result;
}


template<class T>
_DPlane<T> _DPlane<T>::rotate_image(double angle, T bg_val) const
{
  int new_rows = int(ceil(sqrt(_DMatrix<T>::rows()*_DMatrix<T>::rows()+_DMatrix<T>::cols()*_DMatrix<T>::cols())));
  int new_cols = new_rows;

  _DPlane<T> new_I(new_rows, new_cols);
  new_I = T(bg_val);

  int row_half = (_DMatrix<T>::rows()/2);
  int col_half = (_DMatrix<T>::cols()/2);

  int new_row_half = (new_rows/2);
  int new_col_half = (new_cols/2);

  float c = cos(angle);
  float s = sin(angle);


  T *cp = new_I[0];
  for(int i=0; i<new_rows; ++i)
    {
      float i2 = (c*(i-new_row_half)+s*(-1-new_col_half)) + row_half;
      float j2 = (-s*(i-new_row_half)+c*(-1-new_col_half)) + col_half;

      for(int j=0; j<new_cols; ++j, ++cp)
	{
      	  i2 += s; j2 += c;

	  if(i2 < 0 || j2 < 0 || i2>=_DMatrix<T>::rows()-1 || j2>=_DMatrix<T>::cols()-1)
	  	    continue;

	  int low_i = int(i2);
	  int low_j = int(j2);
	  float i_diff = i2-low_i;
	  float j_diff = j2-low_j;
	  
	  T *cp1 = (*this)[low_i];
	  T *cp2 = (*this)[low_i+1];
	  *cp = T(((1-i_diff)*(1-j_diff)*cp1[low_j]) +
		  ((1-i_diff)*j_diff*cp1[low_j+1]) + 
		  (i_diff*(1-j_diff)*cp2[low_j]) +
		  (i_diff * j_diff * cp2[low_j+1]));
	  
	}
    }

  return new_I;
}


template<class T>
_DPlane<T> _DPlane<T>::rotate_image_nn(double angle) const
{
  int new_rows = int(ceil(sqrt(_DMatrix<T>::rows()*_DMatrix<T>::rows()+_DMatrix<T>::cols()*_DMatrix<T>::cols())));
  int new_cols = new_rows;

  _DPlane<T> new_I(new_rows, new_cols);
  new_I = T(1e100);

  int row_half = (_DMatrix<T>::rows()/2);
  int col_half = (_DMatrix<T>::cols()/2);

  int new_row_half = (new_rows/2);
  int new_col_half = (new_cols/2);

  float c = cos(angle);
  float s = sin(angle);


  T *cp = new_I[0];
  for(int i=0; i<new_rows; ++i)
    {
      float i2 = (c*(i-new_row_half)+s*(-1-new_col_half)) + row_half;
      float j2 = (-s*(i-new_row_half)+c*(-1-new_col_half)) + col_half;

      for(int j=0; j<new_cols; ++j, ++cp)
	{
      	  i2 += s; j2 += c;

	  if(i2 < 0 || j2 < 0 || i2>=_DMatrix<T>::rows()-1 || j2>=_DMatrix<T>::cols()-1)
	    continue;

	  int rnd_i = (int) round(i2);
	  int rnd_j = (int) round(j2);

	  *cp = T((*this)[rnd_i][rnd_j]);
	  
	}
    }

  return new_I;
}

template<class T>
_DPlane<T> _DPlane<T>::binary_thin(void) const
{
  static bool lookup_table[256];

  make_thin_lookup_table(lookup_table);

  _DPlane<T> result = (*this);


  for(int i=1; i<_DMatrix<T>::rows()-1; i++)
    for(int j=1; j<_DMatrix<T>::cols()-1; j++)
      {
	int val=0, n=1;

	for(int i2=-1; i2<=1; i2++)
	  for(int j2=-1; j2<=1; j2++)
	    {
	      if(i2 == 0 && j2 == 0)
		continue;

	      if(result[i+i2][j+j2])
		val += n; 

	      n=n*2;
	    }

	if(!lookup_table[val])
	  result[i][j] = 0;
      }

  return result;
}



#define DECLARE(x) \
  template class _DPlane<x>;

DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)



/*
DPlane DPlane::rescale(const DPlane &img, double row_scale, double col_scale) const
{
  DPlane result1(img._DMatrix<T>::rows()*row_scale, img._DMatrix<T>::cols());

  if(row_scale < 1.0)
    {
      for(int i=0; i<img._DMatrix<T>::rows()*row_scale; i++)
	for(int j=0; j<img._DMatrix<T>::cols(); j++)
	  result1[i][j] = img[i/row_scale][j];
    }
  else if(row_scale > 1.0)
    {
      
    }
  else
    result1 = img;

  DPlane result2(img._DMatrix<T>::rows()*row_scale, img._DMatrix<T>::cols()*col_scale);

  if(col_scale < 1.0)
    {
      for(int i=0; i<img._DMatrix<T>::rows(); i++)
	for(int j=0; j<img._DMatrix<T>::cols()*col_scale; j++)
	  result2[i][j] = result1[i][j/col_scale];
    }
  else if(col_scale > 1.0)
    {
    }
  else
    result2 = result1;

  return result2;
}
*/

