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
#ifndef _DMATRIX_H_
#define _DMATRIX_H_

//#define GSL_SUPPORT
#ifdef GSL_SUPPORT
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#endif

#include <iostream>
#include "DPoint.h"
#include "DRect.h"
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <algorithm>

//#define MIN(a,b) ((a<b)?a:b)
//#define MAX(a,b) ((a>b)?a:b)

inline double MIN(double a, double b)
{
  return (a<b)?a:b;
}

inline float MIN(float a, float b)
{
  return (a<b)?a:b;
}

inline int MIN(int a, int b)
{
  return (a<b)?a:b;
}

inline int MAX(int a, int b)
{
  return (a>b)?a:b;
}

inline float MAX(float a, float b)
{
  return (a>b)?a:b;
}

inline double MAX(double a, double b)
{
  return (a>b)?a:b;
}


template<class T>
class _DMatrix
{
    public:
        _DMatrix(int _rows, int _cols);
        _DMatrix(int _rows, int _cols, const T *array);
        _DMatrix();
        _DMatrix(const _DMatrix<T> &other);

        ~_DMatrix();

        _DMatrix<T> operator*(const _DMatrix<T> &other) const;
        _DMatrix<T> operator+(const _DMatrix<T> &other) const;
        _DMatrix<T> operator-(const _DMatrix<T> &other) const;
        _DMatrix<T> &operator=(const _DMatrix<T> &other);
        _DMatrix<T> operator+(T value) const;
        _DMatrix<T> operator-(T value) const;
        _DMatrix<T> operator*(T value) const;
        _DMatrix<T> operator/(T value) const;
	_DMatrix<T> &operator=(T other);

        T *operator[](int row) const;
	T operator[](DPoint &pt) const;

        _DMatrix<T> transpose() const;

        int rows() const { return _rows; }
        int cols() const { return _cols; }
	DPoint size() const { return DPoint(rows(), cols()); }

//	DRect size() const { return DRect(0,0,rows(),cols()); }

	template<class T2>
        friend _DMatrix<T2> operator+(T2 value, const _DMatrix<T2> &other);
	template<class T2>
        friend _DMatrix<T2> operator-(T2 value, const _DMatrix<T2> &other);
	template<class T2>
        friend _DMatrix<T2> operator*(T2 value, const _DMatrix<T2> &other);
	template<class T2>
        friend std::istream &operator>>(std::istream &is, _DMatrix<T2> &matrix);
	template<class T2>
        friend std::ostream &operator<<(std::ostream &is, const _DMatrix<T2> &matrix);	
	template<class T2>
        friend bool same_size(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2);
	template<class T2>
	friend _DMatrix<T2> operator-(const _DMatrix<T2> &m);

	_DMatrix<T> LU_factor(void);
        _DMatrix<T> extract(const DRect &rect) const;
        _DMatrix<T> extract_row(int row) const;
        _DMatrix<T> extract_col(int col) const;

        _DMatrix<T> reshape(int new_rows, int new_cols);

        T length( void ) const
        { 
	  T result = 0; for(int i=0; i<rows(); i++)
	    for(int j=0; j<cols(); j++)
	      result += (*this)[i][j]* (*this)[i][j];
	  
	  return T(sqrt(double(result)));
	}

	template<class T2>
	friend _DMatrix<T2> pointwise_multiply(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2);
	template<class T2>
	friend _DMatrix<T2> pointwise_divide(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2);
	template<class T2>
	friend _DMatrix<T2> sqrt(const _DMatrix<T2> &m);
	template<class T2>
	  friend _DMatrix<T2> fabs(const _DMatrix<T2> &m);
	template<class T2>
	  friend _DMatrix<T2> pointwise_min(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2);
	template<class T2>
	  friend _DMatrix<T2> pointwise_max(const _DMatrix<T2> &m1, const _DMatrix<T2> &m2);

        void set_submatrix(const DPoint &pt, const _DMatrix<T> &in);
        void set_row(int row, const _DMatrix<T> &in);
        void set_col(int col, const _DMatrix<T> &in);

	void swap_rows(int row1, int row2);

	T median() const;
	T mean() const;
	T sum() const;

        _DMatrix<T> operator<(T value);
        _DMatrix<T> operator>(T value);
        _DMatrix<T> operator==(T value);

	_DMatrix<T> covariance(void);
        _DMatrix<T> means(void);


#ifdef GSL_SUPPORT
	_DMatrix<T> inverse();
	T determinant();
	std::pair<_DMatrix<T>, _DMatrix<T> > eigen();

	template<class T2>
	friend gsl_matrix *DMatrix_to_gsl(_DMatrix<T2> &dm);
	template<class T2>
	friend _DMatrix<T2> gsl_to_DMatrix(gsl_matrix *gm);
	template<class T2>
	friend _DMatrix<T2> gsl_to_DMatrix(gsl_vector *gm);
#endif


    protected:
        void deallocate_storage();
        void initialize_storage();

        T **data;
        T *data_area;
        int _rows, _cols;
};

typedef _DMatrix<double> DMatrix;

#define FOR_SCAN(plane, var, type) \
  type *var = plane[0];							\
  for(int ___i___=0, ___max___=plane.rows()*plane.cols(); i<___max___; ++i)



#endif
