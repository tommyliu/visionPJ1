//
//  k-Fan matching and training code, as described in the paper:
//
//    Crandall, Felzenszwalb, Huttenlocher, "Spatial 
//      Priors for Part-based Recognition using Statistical Models," CVPR 2005.
//     
//
//  Author: David Crandall, 2003-2005
//          crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
// 
// 
#include <Appear.h>

using namespace std;

#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)

template<class T>
_DMultiDMatrix<T> Appear<T>::run_model(_DPlane<unsigned char> &edge_map)
{
  _DMultiDMatrix<T> result = run_model(edge_map[0], edge_map.rows(), edge_map.cols());

  return result;
}


template<class T>
_DMultiDMatrix<T> Appear<T>::run_model(unsigned char *edge_map, int rows, int cols)
{
  //  printf("in run_model\n");

  _DMultiDMatrix<T> result(3, part_count, rows, cols);
  result=0;

  T Z[part_count];
  for(int p=0; p<part_count; p++)
    {
      Z[p]=0;
      T *Z_cp = Z+p;
      T *log_p_cp = (log_fg_model.get(p,0))[0];
      T log_bg_model_0_0 = log_bg_model[0][0];

      int f_size = log_fg_model.rows() * log_fg_model.cols();

      for(int i=0; i<f_size; i++, log_p_cp++)
	(*Z_cp) += (*log_p_cp) - log_bg_model_0_0;
    }

  for(int p=0; p<part_count; p++)
    {
      int f_rows = log_fg_model.rows(), f_rows2 = f_rows/2;
      int f_cols = log_fg_model.rows(), f_cols2 = f_cols/2;
      int f_size = f_rows * f_cols;
      _DMatrix<T> &result_p = result.get(p);

      for(int i=f_rows2; i<rows-f_rows2; i++)
	{
	  T *cp = result_p[i]+f_cols2;

	  for(int j=f_cols2; j<cols-f_cols2; j++)
	    *(cp++) = Z[p];
	}

      unsigned char *edge_map_cp = edge_map;

      for(int i=0; i<rows; i++)
	{
	  for(int j=0; j<cols; j++)
	    {
	      int this_dir = *(edge_map_cp++);
	      if(!this_dir)
		continue;

	      int high_bound = MIN(f_rows2,rows-i-1-f_rows2);
	      int low_bound = MAX(-f_rows2,-i-1+f_rows2 );
	      int high_bound1 = MIN(f_cols2,cols-j-1-f_cols2);
	      int low_bound1 = MAX(-f_cols2,-j-1+f_cols2);
	      
	      T *thislogp = (log_back_probs_diff.get(p,this_dir))[0];

	      for(int k=high_bound; k>low_bound; k--)
		{
		  T *out_cp = result_p[i+k]+(j+high_bound1);
		  T *thislogp_cp = thislogp + (f_rows2-k)*f_cols + 
		    f_cols2 - high_bound1;

		  for(int l=high_bound1; l>low_bound1; l--, out_cp--, thislogp_cp++)
		    {
		      *out_cp += *thislogp_cp;
		    }
		}
	    }
	}
    }
    
  return result;
}


template<class T>
unsigned char *AppearModel<T>::embed_img(unsigned char *in_im, int rows, int cols,
				 int new_rows, int new_cols)
{
  unsigned char *new_img = new unsigned char[new_rows*new_cols];

  int r_border = (new_rows-rows)/2;
  int c_border = (new_cols-cols)/2;

  memset(new_img, 0, new_rows*new_cols);


  for(int i=0; i<rows; i++)
    {
      unsigned char *ni_cp = &(new_img[(r_border+i)*new_cols + c_border]);
      unsigned char *cp = in_im+i;

      for(int j=0; j<cols; j++, ni_cp++, cp+=rows)
	*ni_cp = *cp;
    }
      
  return new_img;

}             
 
/*
DMultiDMatrix Appear::run_model_rotation_invariant_exhaustive(unsigned char *img, int rows, int cols, double min_angle, double max_angle,
							      double angle_step)
{
  DPlane in_img(rows, cols);

  unsigned char *in_cp=img;
  double *out_cp=in_img[0];
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++, out_cp++, in_cp++)
      {
	*out_cp = *in_cp;
      }

  return run_model_rotation_invariant_exhaustive(in_img, min_angle, max_angle, angle_step);
}


template<class T>
_DMultiDMatrix<T> Appear<T>::run_model_rotation_invariant_exhaustive(const DPlane &img, T min_angle, T max_angle, T angle_step)
{
  int angle_no=0;

  _DMultiDMatrix<T> best_sofar(3, part_count, edge_map.rows(), edge_map.cols()); 
  best_sofar = 1e100;

  for(T angle = min_angle; angle <= max_angle; angle += angle_step)
    {
      _DPlane<T> in_plane(edge_map);
      in_plane = in_plane.rotate_image_nn(-angle);
 
      _DMultiDMatrix<T> img = run_model(in_plane);
      int height = in_plane.rows(), width = in_plane.cols();

      for(int p=0; p<part_count; p++)
	{
	  cout << p << " --- " << angle << endl;
	  _DPlane<T> plane(img.get(p));
	  _DPlane<T> I4 = plane.rotate_image(angle);
	  
	  int r1 = edge_map.rows(), c1 = edge_map.cols();
	  int r2 = I4.rows(), c2 = I4.cols();
	  
	  int new_rows_half = height/2;
	  int new_cols_half = width/2;
	  
	  int rr1 = ((r2-r1)/2);
	  int cc1 = ((c2-c1)/2);
	  
	  int rr2 = ((height-r1)/2);
	  int cc2 = ((width-c1)/2);
      
	  _DPlane<T> im = I4.extract(DRect(rr1, cc1, rr1+r1-1, cc1+c1-1));

	  best_sofar.get(p) = pointwise_min(im, best_sofar.get(p));
	}

      angle_no++;
    }

  return best_sofar;
}
*/


/*
DMultiDMatrix Appear::run_model_rotation_invariant(const DPlane &img)
{
  DPlane H(img.rows(), img.cols()), V(img.rows(), img.cols());

  for(int i=1; i<img.rows()-1; i++)
    for(int j=1; j<img.cols()-1; j++)
      {
	H[i][j] = img[i+1][j+1] - img[i-1][j-1];
	V[i][j] = img[i+1][j-1] - img[i-1][j+1];
      }

  for(int i=1; i<img.rows()-1; i++)
    for(int j=1; j<img.cols()-1; j++)

	-(1/2*c-1/2*a-1/2*(c^2-2*a*c+a^2+4*b^2)^(1/2))/b
      }


}
*/


template<class T>
_DMultiDMatrix<T> Appear<T>::run_model_with_rotation(unsigned char *img, int rows, int cols, double min_angle, double max_angle,
					      double angle_step)
{
  _DPlane<unsigned char> in_img(rows, cols);

  unsigned char *in_cp=img;
  unsigned char *out_cp=in_img[0];
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++, out_cp++, in_cp++)
      {
	*out_cp = *in_cp;
      }

  return run_model_with_rotation(in_img, min_angle, max_angle, angle_step);
}

template<class T>
_DMultiDMatrix<T> Appear<T>::run_model_with_rotation(_DPlane<unsigned char> &edge_map, double min_angle, double max_angle,
					      double angle_step)
{
  cout << "in run model rotat" << endl;
  int rot_count = int((max_angle-min_angle) / angle_step)+1;
  _DMultiDMatrix<T> cm(4, rot_count, part_count, edge_map.rows(), edge_map.cols()); 

  int angle_no=0;
  for(double angle = min_angle; angle <= max_angle; angle += angle_step)
    {
      _DPlane<unsigned char> in_plane(edge_map);
      in_plane = in_plane.rotate_image_nn(-angle);
 
      _DMultiDMatrix<T> img = run_model(in_plane);
      int height = in_plane.rows(), width = in_plane.cols();

      for(int p=0; p<part_count; p++)
	{
	  cout << p << " --- " << angle << endl;
	  _DPlane<T> plane(img.get(p));
	  _DPlane<T> I4 = plane.rotate_image(angle);
	  
	  int r1 = edge_map.rows(), c1 = edge_map.cols();
	  int r2 = I4.rows(), c2 = I4.cols();
	  
	  int new_rows_half = height/2;
	  int new_cols_half = width/2;
	  
	  int rr1 = ((r2-r1)/2);
	  int cc1 = ((c2-c1)/2);
	  
	  int rr2 = ((height-r1)/2);
	  int cc2 = ((width-c1)/2);
      
	  _DPlane<T> im = I4.extract(DRect(rr1, cc1, rr1+r1-1, cc1+c1-1));

	  cm.get(angle_no, p) = im;
	}
      angle_no++;
    }
  
  return cm;
}


#define DECLARE(x) \
  template class Appear<x>; \
  template class AppearModel<x>;



DECLARE(double)
DECLARE(float)
