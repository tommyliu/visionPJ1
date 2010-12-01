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
#ifndef __DMULTIDMATRIX_H_
#define __DMULTIDMATRIX_H_

#include "DRect.h"
#include "DMatrix.h"

template<class T>
class _DMultiDMatrix
{
 public:
  _DMultiDMatrix() 
    {
      dimensions = 0;
      dim_sizes = 0;
      matrices = 0;
      mat_count=0;
    }

  _DMultiDMatrix(int _dim_count, int p1, int p2, int p3)
    {
      assert(_dim_count == 3);
	
      int dims[3];
      dims[0] = p1;
      dims[1] = p2;
      dims[2] = p3;

      allocate_storage(3, dims);
    }
  
  _DMultiDMatrix(int _dim_count, int p1, int p2, int p3, int p4)
    {
      assert(_dim_count == 4);
	
      int dims[4];
      dims[0] = p1;
      dims[1] = p2;
      dims[2] = p3;
      dims[3] = p4;

      allocate_storage(4, dims);
    }

  _DMultiDMatrix(int _dim_count, int *_dims)
    {
            assert(_dim_count > 2);

      allocate_storage(_dim_count, _dims);
    }
  

  void allocate_storage(int _dim_count, int *_dims)
  {
      dimensions = _dim_count;
      dim_sizes = new int[dimensions];
      for(int i=0; i<dimensions; i++)
	dim_sizes[i] = _dims[i];

      mat_count = 1;
      for(int i=0; i<dimensions-2; i++)
	mat_count *= _dims[i];

      matrices = new _DMatrix<T> *[mat_count];

      for(int i=0; i<mat_count; i++)
	matrices[i] = new _DMatrix<T>(_dims[dimensions-2], _dims[dimensions-1]);    
  }

  
  ~_DMultiDMatrix()
    {
      if(dim_sizes)
	delete[] dim_sizes;
      
      if(matrices)
	{
	  for(int i=0; i<mat_count; i++)
	    delete matrices[i];
	  
	  delete[] matrices;
	}
    }
  
    _DMultiDMatrix(const _DMultiDMatrix<T> &other)
      {
	matrices = 0;
	dim_sizes = 0;

	*this = other;
      }

    //    DMultiDMatrix &operator[](int plane) const
    //    {
    //       return *next_level[plane];
    //    }

    _DMultiDMatrix<T> &operator=(const T D)
      {
	for(int i=0; i<mat_count; i++)
	  get(i) = D;
      }

    _DMultiDMatrix<T> &operator=(const _DMultiDMatrix<T> &other)
      {
       if(matrices)
	 {
	   for(int i=0; i<mat_count; i++)
	     delete matrices[i];

	   delete[] matrices;
	 }

       if(dim_sizes)
	 delete[] dim_sizes;

       mat_count = other.planes();
       matrices = new _DMatrix<T> *[mat_count];
       dimensions = other.get_dimensions_count();

       for(int i=0; i<mat_count; i++)
	 matrices[i] = new _DMatrix<T>(other.get(i));

       dim_sizes = new int[dimensions];
       for(int i=0; i<dimensions; i++)
	 dim_sizes[i] = other.get_dims(i);

      }

    int get_dimensions_count() const
    {
      return dimensions;
    }

    _DMatrix<T> &get(int c1) const
    {
      return *(matrices[c1]);
    }

    _DMatrix<T> &get(int c1, int c2) const
    {
      return *(matrices[c1*dim_sizes[1]+c2]);
    }

    template<class T2>
      friend _DMultiDMatrix<T2> pointwise_min(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2);
    template<class T2>
      friend _DMultiDMatrix<T2> pointwise_max(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2);


    int get_dims(int i) const
    {
      return dim_sizes[i];
    }

    template<class T2>
      friend bool same_size(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2);


    int rows() const { return matrices[0]->rows(); }
    int cols() const { return matrices[0]->cols(); }
    int planes() const { return mat_count; }

    friend std::istream &operator>>(std::istream &is, _DMultiDMatrix<T> &matrix)
    {
      int dimensions;
      is >> dimensions;

      if(dimensions)
	{
	  int dims[dimensions];
	  for(int i=0; i<dimensions; i++)
	    is >> dims[i];
	  
	  matrix = _DMultiDMatrix<T>(dimensions, dims);
	  
	  for(int i=0; i<matrix.planes(); i++)
	    is >> matrix.get(i);
	}
      return is;
    }

    friend std::ostream &operator<<(std::ostream &os, _DMultiDMatrix<T> &matrix)
    {
      os << matrix.dimensions << std::endl;

      for(int i=0; i<matrix.dimensions; i++)
	os << matrix.dim_sizes[i] << " ";

      for(int i=0; i<matrix.planes(); i++)
	os << matrix.get(i);

      return os;
    }

    

  protected:
    int dimensions;
    int mat_count;

    int *dim_sizes;
    _DMatrix<T> **matrices;
};


typedef _DMultiDMatrix<double> DMultiDMatrix;

template class _DMultiDMatrix<double>;
template class _DMultiDMatrix<float>;
template class _DMultiDMatrix<int>;
template class _DMultiDMatrix<short>;
template class _DMultiDMatrix<char>;


#endif

