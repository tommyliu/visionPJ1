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
#ifndef _DINTEGRAL_IMAGE_H_
#define _DINTEGRAL_IMAGE_H_

#include "DPlane.h"
#include "DRect.h"

class DIntegralImage : public DPlane
{
 public:
  DIntegralImage(DPlane &plane)
  {
    construct_integral_image(plane);
  }

  double get_sum(DRect &rect);

 protected:
  void construct_integral_image(DPlane &plane);
  DPlane integral_image;
};

#endif
