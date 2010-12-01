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
#include <vector>
#include "DPlane.h"
#include "DImage.h"

template<class T, class T2>
  class _DistanceTransform
{
 public:
  _DistanceTransform() {}

  std::pair<_DPlane<T>,std::pair<_DPlane<T2>, _DPlane<T2> > > do_transform_2d(const _DPlane<T> &src_img, T scaling=1);
  std::pair<_DPlane<T>,std::pair<_DPlane<T2>, _DPlane<T2> > > do_transform_2d(const _DPlane<T> &src_img, T x_scaling, T y_scaling);

  // distance transform, with arbitrary (i.e. possibly non-diagonal)
  // covariance matrix
  std::pair<_DPlane<T>,std::pair<_DPlane<T2>, _DPlane<T2> > > do_transform_2d(const _DPlane<T> & in_im, _DMatrix<T> &sigma);


  DImage do_transform_3d(const DImage & in_im, _DMatrix<T> &sigma, 
					    T s_z);
};

typedef _DistanceTransform<double, double> DistanceTransform;
