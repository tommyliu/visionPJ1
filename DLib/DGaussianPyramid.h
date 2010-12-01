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
#ifndef __GAUSSIANPYR_H_
#define __GAUSSIANPYR_H_

#include "DRect.h"
#include "DPlane.h"
#include "DIntegralImage.h"
#include "DGaussian.h"

class DGaussianPyramid 
{
  public:
    DGaussianPyramid(int _scales, const DPlane &in_plane);
    ~DGaussianPyramid();

    DGaussianPyramid(const DGaussianPyramid &other);

    DPlane &operator[](int plane) const
    {
       return plane_data[plane];
    }

    DGaussianPyramid &operator=(const DGaussianPyramid &other);

    int scales() const { return scale_count; }

    void compute_gaussian_pyramid(const DPlane &in_image) const;

  protected:
    int scale_count;
    DPlane *plane_data;

};

#endif
