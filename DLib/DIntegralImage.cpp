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
#include "DIntegralImage.h"


void DIntegralImage::construct_integral_image(DPlane &plane)
{
  integral_image = DPlane(plane.rows(), plane.cols());

  // horizontal pass
  double *in_ptr = plane[0];
  double *ii_ptr = integral_image[0];
  for(int i=0; i<plane.rows(); i++)
    {
      double sum_so_far = 0;
      for(int j=0; j<plane.cols(); j++, in_ptr++, ii_ptr++)
	sum_so_far = *ii_ptr = sum_so_far + *in_ptr;
    } 

  // vertical pass
  int col_count = plane.cols();
  for(int j=0; j<plane.cols(); j++)
    {
      ii_ptr = integral_image[0]+j;
      double sum_so_far = 0;
      for(int i=0; i<plane.rows(); i++, ii_ptr += col_count)
	sum_so_far = *ii_ptr += sum_so_far;
    }

}


double DIntegralImage::get_sum(DRect &rect)
{
  double result;
  double top_region, left_region, corner_region;

  if(rect.top() == 0)
    top_region = 0;
  else
    top_region = integral_image[rect.top()-1][rect.right()];

  if(rect.left() == 0)
    left_region = 0;
  else
    left_region = integral_image[rect.bottom()][rect.left()-1];

  if(rect.top() == 0 || rect.left() == 0)
    corner_region = 0;
  else
    corner_region = integral_image[rect.top()-1][rect.left()-1]; 

  result = integral_image[rect.bottom()][rect.right()] +
            corner_region - top_region - left_region;

//  result = integral_image[rect.bottom_right()] +
//    2.0 * integral_image[rect.top_left()] -
//    integral_image[rect.top_right()] -
// integral_image[rect.bottom_left()];

  return result;
}
