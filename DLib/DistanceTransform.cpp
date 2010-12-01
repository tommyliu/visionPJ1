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
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "DistanceTransform.h"
using namespace std;

/* Implementation of fast squared Euclidean distance transform algorithm
   using amortized algorithm for lower envelope of quadratics.  

   For a description see
   www.cs.cornell.edu/~dph/matchalgs/iccv2003-tutorial.pdf */

// hacked for C++ by crandall, 9/2003

template<class T,class T2>
pair<_DPlane<T>,pair<_DPlane<T2>, _DPlane<T2> > > 
_DistanceTransform<T,T2>::do_transform_2d
(const _DPlane<T> & in_im, T scaling) 
{
  return do_transform_2d(in_im, scaling, scaling);
}

template<class T,class T2>
pair<_DPlane<T>,pair<_DPlane<T2>, _DPlane<T2> > > _DistanceTransform<T,T2>::do_transform_2d(const _DPlane<T> & in_im, T _scaling_x, T _scaling_y) 
{
  float scaling_x = _scaling_x;
  float scaling_y = _scaling_y;

  _DPlane<T> im = in_im;

  // records closest pixel
  _DPlane<T2> closest_rows(in_im.rows(), in_im.cols());
  _DPlane<T2> closest_cols(in_im.rows(), in_im.cols());

  int width = im.cols();
  int height = im.rows();
  int *z, *v;
  float *vref;
  int k;
  float s;
  float sp;
  int x, y;

  z = (int *)malloc(sizeof(int)*(max(width, height)+1));
  v = (int *)malloc(sizeof(int)*(max(width, height)));
  vref = (float *)malloc(sizeof(float)*(max(width, height)));
  if ((z == NULL) || (v == NULL) || (vref == NULL)) {
    assert(0);
  }


  static int first=1;
  static float lut1[2000];
  //  if(first)
    {
  float sx_inv = 1.0/scaling_x;
      for(int i=-height,j=0; i<height; i++, j++) //y
	  lut1[j] = 1.0/(2.0 * i) * sx_inv;
      first=0;
    }

    static int _first=1;
    static float *lut2[2000];
    if(_first)
      {
	for(int i=0; i<height; i++)
	{
	  lut2[i] = new float[height];

	for(int j=0; j<height; j++)
	  {
	    lut2[i][j] = i*i-j*j;
	  }
	}
	_first=0;
      }


  /* do x transform */
  for (y = 0; y < height; y++) {
    k = 0;  /* Number of boundaries between parabolas */
    z[0] = 0;   /* Indexes of locations of boundaries,
		   order by increasing x */
    z[1] = width;    
    v[0] = 0;     /* Indexes of locations of visible parabola bases,
		     ordered by increasing x */

    T *im_y=im[y];
    T2 *closest_rows_y = closest_rows[y];
    T2 *closest_cols_y = closest_cols[y];

    int v_k = v[k], v_k_2 = v[k]*v[k];

    for (x = 1; x < width; x++) {
      // # of times here = 1.0
      do {
	/* compute Vornoi border: intersection of parabola at x
	   with rightmost currently visible parabola */
	//	s = ((im[y][x] + x*x) - (im[y][v[k]] + v[k]*v[k])) /
	//	  (2 * (x - v[k]));
	//	s = ((im_y[x] + scaling_x*x*x) - (im_y[v[k]] + scaling_x*v[k]*v[k])) /
	//	  (2 * scaling_x * (x - v[k]));


	// count 1.1
	s = (im_y[x] - im_y[v_k] + scaling_x*(x*x - v_k_2)) *
	  lut1[x-v_k+height]; // * sx_inv;
	  // / (2 * scaling_x * (x - v_k));
	  
	sp = ceil(s); // floor(ceil(s));

	/* case one: intersection is to the right of the array, 
	   so this parabola is not visible (and nothing to do) */
	if (sp >= width)
	  { 
	    // count 0.8
	  break;
	  }

	/* case two: intersection is at larger x than rightmost current
	       intersection, so this parabola is visible on the right (add
	       it to the end) */
	if (sp > z[k]) {
	  z[k+1] = int(sp);
	  z[k+2] = width;
	  v[k+1] = x;
	  k++; 
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	  // count 0.2
	  break;
	}

	/* case three: intersection is at smaller x than the
	     rightmost current intersection, so this parabola hides the
	       rightmost parabola (remove the rightmost parabola, if there
	         are still remaining potentially visible parabolas iterate to
		 complete the addition of this parabola). */

	if (k == 0) {
	  v[k] = x;
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	  // count 0.002
	  break;
	} else {
	  z[k] = width;
	  k--;
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	  // count 0.09
	}
      } while (1);
      
    }

    /* compute transform values from visible parabolas */
    
    /* get value of input image at each parabola base */
    for (x = 0; x <= k; x++) {
      vref[x] = im_y[v[x]];
    }
    k = 0;

    /* iterate over pixels, calculating value for closest parabola */
    v_k=v[k];
    for (x = 0; x < width; x++) {
      if (x == z[k+1])
	k++, v_k=v[k];
      im_y[x] = vref[k] + (v_k-x)*(v_k-x) * scaling_x;

      closest_rows_y[x] = y;
      closest_cols_y[x] = v_k;
    }

  }


  float sy_inv = 1.0/scaling_y;
  for(int i=-height,j=0; i<height; i++, j++) //y
    lut1[j] = 1.0/(2.0 * i) * sy_inv;
  first=0;  
  
  


  _DPlane<T2> closest_rows2(closest_rows);
  _DPlane<T2> closest_cols2(closest_cols);



  /* do y transform - analogous computation in y-direction applied to
     result of x-transform */
  for (x = 0; x < width; x++) {
    k = 0;
    z[0] = 0;
    z[1] = height;
    v[0] = 0;

    int v_k = v[k];
    int v_k_2 = v[k]*v[k];
    T *im_x_cp=im[1]+x;

    for (y = 1; y < height; y++, im_x_cp+=width) {

      do {
	/* compute vornoi border */
	float s1 = *im_x_cp;
	s1+=(y*y-v_k_2)*scaling_y;
	s1 -= im[v_k][x];
	//	float s2 = (2 * (y - v_k) * scaling_y);
	  
	float s=s1*lut1[y-v_k+height]; //*sy_inv;
	sp = ceil(s);

	/* case one */
	if (sp >= height)
	  break;

	/* case two */
	if (sp > z[k]) {
	  z[k+1] = int(sp);
	  z[k+2] = height;
	  v[k+1] = y;
	  k++;
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	  break;
	}

	/* case three */
	if (k == 0) {
	  v[0] = y;
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	  break;
	} else {
	  z[k] = height;
	  k--;
	  v_k = v[k];
	  v_k_2 = v_k*v_k;
	}
      } while (1);
      
    }
    
    for (y = 0; y <= k; y++) {
      vref[y] = im[v[y]][x];
    }
    k = 0;

    T2 *closest_cols2__x = closest_cols2[0]+x;
    T2 *closest_rows2__x = closest_rows2[0]+x;
    v_k = v[k];
    T2 *closest_rows_in__x = closest_rows[v_k]+x;
    T2 *closest_cols_in__x = closest_cols[v_k]+x;
    //    double *im_y = im[0]+x;

    for (y = 0; y < height; y++, closest_cols2__x+=width, 
	   closest_rows2__x+=width) //, im_y += width) 
      {
	if (y == z[k+1])
	  {
	    k++;
	    closest_cols_in__x += width*(v[k]-v_k);
	    closest_rows_in__x += width*(v[k]-v_k);
	    v_k = v[k];
	  }

	im[y][x] = vref[k] + (v_k-y)*(v_k-y)*scaling_y;
	
	*closest_cols2__x = *closest_cols_in__x;
	*closest_rows2__x = *closest_rows_in__x;
      }
    
  }


  
  free(z);
  free(v);
  free(vref);
  return pair<_DPlane<T>, pair<_DPlane<T2>, _DPlane<T2> > >
    (im,pair<_DPlane<T2>, _DPlane<T2> >(closest_rows2,closest_cols2));
}

  
// distance transform, with arbitrary (i.e. possibly non-diagonal)
// covariance matrix
template<class T,class T2>
pair<_DPlane<T>,pair<_DPlane<T2>, _DPlane<T2> > > _DistanceTransform<T,T2>::do_transform_2d(const _DPlane<T> & orig_im, _DMatrix<T> &sigma)
{
  int rotate;
  double angle=0, scaling_x, scaling_y;
  _DPlane<T> in_im = orig_im;
  
  if(sigma[0][1] == 0 && sigma[1][0] == 0)
    rotate=0;
  else
    rotate=1;


  if(rotate)
    {
#ifdef GSL_SUPPORT
      pair<_DMatrix<T>, _DMatrix<T> > res = sigma.eigen();

      if(res.second[0][1] ==0)
	angle=0;
      else
	angle = -atan(res.second[1][1]/res.second[0][1]);
    
      scaling_x = 1/res.first[0][1];
      scaling_y = 1/res.first[0][0];
  
      in_im = in_im.rotate_image(-angle);
#else
      assert(0);
#endif

      //     ind = find(isinf(in_im));
      // in_im(ind) = 10e100;
    }
  else
    {
      scaling_x = 1/sigma[0][0];
      scaling_y = 1/sigma[1][1];
    }

  int height = in_im.rows();
  int width = in_im.cols();
  //  cout << angle << endl;
  //  cout << in_im;

  pair<_DPlane<T>, pair<_DPlane<T2>, _DPlane<T2> > > result = 
    do_transform_2d(in_im, scaling_x, scaling_y);
 
  _DPlane<T> &dt_result = result.first;
  _DPlane<T2> &closest_x = result.second.second;
  _DPlane<T2> &closest_y = result.second.first;
 
  if(rotate)
    {
      _DPlane<T> I4 = dt_result.rotate_image(angle);

      int r1 = orig_im.rows(), c1 = orig_im.cols();
      int r2 = I4.rows(), c2 = I4.cols();

      int new_rows_half = height/2;
      int new_cols_half = width/2;

      int rr1 = ((r2-r1)/2);
      int cc1 = ((c2-c1)/2);

      int rr2 = ((height-r1)/2);
      int cc2 = ((width-c1)/2);

      _DPlane<T> im = I4.extract(DRect(rr1, cc1, rr1+r1-1, cc1+c1-1));

      DMatrix R(2,2);
      R[0][0] = R[1][1] = cos(-angle);
      R[0][1] = sin(-angle);
      R[1][0] = -sin(-angle);

      // fix closest_rows, closest_cols here (to compensate for rotation)
      closest_x = closest_x.rotate_image(angle);
      closest_y = closest_y.rotate_image(angle);

      closest_x = closest_x.extract(DRect(rr1, cc1, r1+rr1-1, c1+cc1-1));
      closest_y = closest_y.extract(DRect(rr1, cc1, r1+rr1-1, c1+cc1-1));

      double cos_ = cos(-angle);
      double sin_ = sin(-angle);

      DMatrix pt(2,1);
      for(int i=0; i<r1; i++)
	{
	  T2* closest_y_i = closest_y[i];
	  T2* closest_x_i = closest_x[i];

	  for(int j=0; j<c1; j++)
	    {
	      //	      pt[0][0] = closest_y_i[j]-new_rows_half;
	      //	      pt[1][0] = closest_x_i[j]-new_cols_half;
	    
	      //	      DMatrix new_pt = R * pt;

	      T2 old_row = closest_y_i[j]-new_rows_half;
	      T2 old_col = closest_x_i[j]-new_cols_half;

	      double new_row = cos_ * old_row + sin_ * old_col;
	      double new_col = -sin_ * old_row + cos_ * old_col;

	      closest_y_i[j] = T2(rint(new_row + new_rows_half - rr2));
	      closest_x_i[j] = T2(rint(new_col + new_cols_half - cc2));
	    }
	}

      dt_result = im;

    }

  return(pair<_DPlane<T>, pair<_DPlane<T2>,_DPlane<T2> > >
	 (dt_result, pair<_DPlane<T2>,_DPlane<T2> >(closest_x, closest_y)));

}




#ifdef COMMENTOUT
// 3-d distance transform, with arbitrary (i.e. possibly non-diagonal)
// covariance matrix
template<class T, class T2>
DImage _DistanceTransform<T,T2>::do_transform_3d(const DMultiDMatrix & in_im, _DMatrix<T> &sigma, 
					  T s_z)
{
  int rows = in_im.rows();
  int cols = in_im.cols();
  int planes = in_im.planes();

  pair<DPlane, pair<DPlane, DPlane> > *result = 
    new pair<DPlane, pair<DPlane, DPlane> >[planes];
  DistanceTransform dt;

  for(int k=0; k<planes; k++)
    {
      result[k] = dt.do_transform_2d(in_im.get(k), sigma);
    }

  int width = in_im.cols();
  int height = in_im.rows();
  int *z, *v;
  float *vref;
  int k;
  float s;
  float sp;
  int x, y;

  z = (int *)malloc(sizeof(int)*planes);
  v = (int *)malloc(sizeof(int)*planes);
  vref = (float *)malloc(sizeof(float)*planes);
  if ((z == NULL) || (v == NULL) || (vref == NULL)) {
    assert(0);
  }

  /* do z transform */
  for (int p = 0; p < planes; p++) {
    k = 0;  /* Number of boundaries between parabolas */
    z[0] = 0;   /* Indexes of locations of boundaries,
		   order by increasing x */
    z[1] = width;    
    v[0] = 0;     /* Indexes of locations of visible parabola bases,
		     ordered by increasing x */

    T *im_p=result.get(p);
    T2 *closest_rows_y = closest_rows.get(p);
    T2 *closest_cols_y = closest_cols.get(p);

    for(int i=0; i<height; i++)
      for(int j=0; j<width; j++)
	{
	  do {
	    /* compute Vornoi border: intersection of parabola at x
	       with rightmost currently visible parabola */
	    s = ((im_p[i][j] + scaling_z*p*p) - (im_y[v[k]] + scaling_z*v[k]*v[k])) /
	      (2 * scaling_z * (p - v[k]));
	
	    sp = ceil(s); // floor(ceil(s));
	
	    /* case one: intersection is to the right of the array, 
	       so this parabola is not visible (and nothing to do) */
	    if (sp >= width)
	      { 
		break;
	      }

	    /* case two: intersection is at larger x than rightmost current
	       intersection, so this parabola is visible on the right (add
	       it to the end) */
	    if (sp > z[k]) {
	      z[k+1] = int(sp);
	      z[k+2] = width;
	      v[k+1] = x;
	      k++; 
	      break;
	    }
	    
	    /* case three: intersection is at smaller x than the
	       rightmost current intersection, so this parabola hides the
	       rightmost parabola (remove the rightmost parabola, if there
	       are still remaining potentially visible parabolas iterate to
	       complete the addition of this parabola). */
	    
	    if (k == 0) {
	      v[k] = x;
	      break;
	    } else {
	      z[k] = width;
	      k--;
	    }
	  } while (1);
	  
	}
    
    /* compute transform values from visible parabolas */
    
    /* get value of input image at each parabola base */
    for (x = 0; x <= k; x++) {
      vref[x] = im_y[v[x]];
    }
    k = 0;
    
    /* iterate over pixels, calculating value for closest parabola */
    v_k=v[k];
    for (x = 0; x < width; x++) {
      if (x == z[k+1])
	k++, v_k=v[k];
      im_y[x] = vref[k] + (v_k-x)*(v_k-x) * scaling_x;

      closest_rows_y[x] = y;
      closest_cols_y[x] = v_k;
    }

  }








  DMultiDMatrix result2(3, planes, rows, cols);

  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      {
	for(int k=0; k<planes; k++)
	  {

	    // k=0
	    result2[0][i][j] = min(
				   min(result[0].first[i][j],
				       result[1].first[i][j] + s_z),
				   min(result[2].first[i][j] + 2*s_z,
				       result[3].first[i][j] + s_z));
	    
	    // k=1
	    result2[1][i][j] = min(
				   min(result[0].first[i][j] + s_z,
				       result[1].first[i][j]),
				   min(result[2].first[i][j] + s_z,
				       result[3].first[i][j] + 2*s_z));
	    
	    // k=2
	    result2[2][i][j] = min(
				   min(result[0].first[i][j] + 2*s_z,
				       result[1].first[i][j] + s_z),
				   min(result[2].first[i][j],
				       result[3].first[i][j] + s_z));
	    
	    // k=3
	    result2[3][i][j] = min(
				   min(result[0].first[i][j] + s_z,
				       result[1].first[i][j] + 2*s_z),
				   min(result[2].first[i][j] + s_z,
				       result[3].first[i][j]));
	    
	  }
      }

  delete[] result;  

  return result2;
}
#endif

#define DECLARE(x,y)			\
  template class _DistanceTransform<x,y>; 

DECLARE(double,double);
DECLARE(double, short);
DECLARE(float, short);
