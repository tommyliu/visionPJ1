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
#include <DGaussianPyramid.h>

DGaussianPyramid::DGaussianPyramid(int _scales, const DPlane &in_plane)
{
  plane_data = new DPlane[_scales];
  scale_count = _scales;
  
  for(int i=0; i<_scales; i++)
    plane_data[i] = DPlane(100,100);

  compute_gaussian_pyramid(in_plane);
}

DGaussianPyramid::~DGaussianPyramid()
{
  delete[] plane_data;
}

DGaussianPyramid::DGaussianPyramid(const DGaussianPyramid &other)
{
  *this = other;
}

DGaussianPyramid &DGaussianPyramid::operator=(const DGaussianPyramid &other)
{
  scale_count = other.scales();
  plane_data = new DPlane[other.scales()];
  
  for(int i=0; i<scale_count; i++)
    plane_data[i] = DPlane(other[i]);
  
  return *this;
}

void DGaussianPyramid::compute_gaussian_pyramid(const DPlane &in_image) const
{
  int rows = in_image.rows();
  int cols = in_image.cols();

  DGaussian g_filter(7, 7, 1.0);

  // first plane (original image)
  plane_data[0] = in_image;

  for(int i=1; i<scales(); i++)
    {
      printf("scale %d\n", i);
      DPlane smoothed = plane_data[i-1].convolve(g_filter);

      plane_data[i] = smoothed.subsample(2, 2);
    }
}



