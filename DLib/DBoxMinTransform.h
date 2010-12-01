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
#ifndef __DBOXMINTRANSFORM_H_
#define __DBOXMINTRANSFORM_H_

#include <DImage.h>
#include <DPlane.h>

template<class T>
class DBoxMinTransform
{
 public:
  DBoxMinTransform(const _DMatrix<T> &img, int max_r, int max_c);
  ~DBoxMinTransform();

  void del_prime_d(T *, int, int, int, T **);
  void del_prime_d2(T *, int, int, int, int,  T **);
  _DMatrix<T> do_transform(int d_row, int d_col);
			 
					      

 protected:
  _DMatrix<T> Delta1;
  int max_log_d;
  T ***del_primes_2n;

};


#endif
