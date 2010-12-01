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
#include<DMultiDMatrix.h>

template<class T2>
_DMultiDMatrix<T2> pointwise_min(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2)
{
  assert(same_size(m1, m2));
  
  _DMultiDMatrix<T2> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = pointwise_min(m1.get(i), m2.get(i));
  
  return result;
}

template<class T2>
_DMultiDMatrix<T2> pointwise_max(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2)
{
  assert(same_size(m1, m2));
  
  _DMultiDMatrix<T2> result(m1);
  for(int i=0; i<m1.planes(); i++)
    result.get(i) = pointwise_max(m1.get(i), m2.get(i));
  
  return result;
}

template<class T2>
bool same_size(const _DMultiDMatrix<T2> &m1, const _DMultiDMatrix<T2> &m2)
{
  return m1.rows() == m2.rows() && m1.cols() == m2.cols() && m1.planes() == m2.planes();
}



#define DECLARE(x) \
  template class _DMatrix<x>; \
  template _DMultiDMatrix<x> pointwise_min(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2); \
  template _DMultiDMatrix<x> pointwise_max(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2); \
  template bool same_size(const _DMultiDMatrix<x> &m1, const _DMultiDMatrix<x> &m2);



DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)
