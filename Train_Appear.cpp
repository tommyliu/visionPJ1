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
#include <Train_Appear.h>
using namespace std;

template<class T>
_DMatrix<char> Train_Appear<T>::read_edge_map(const char *fname)
{
  int rows, cols;

  FILE *fp = fopen(fname, "rb");
  fread(&rows, 1, sizeof(int), fp);
  fread(&cols, 1, sizeof(int), fp);

  _DMatrix<char> result(cols, rows);
  
  fread(result[0], rows*cols, sizeof(char), fp);
  fclose(fp);

  return result.transpose();
}


template<class T>
_DMatrix<char> Train_Appear<T>::get_image_patch(const _DMatrix<char> &Q, 
			int row, int col, int tile_rows, int tile_cols)
{
  _DMatrix<char> result(tile_rows, tile_cols);

  int tile_rows_half = tile_rows/2;
  int tile_cols_half = tile_cols/2;

  for(int i=row-tile_rows_half, i2=0; i<row+tile_rows_half; i++, i2++)
    for(int j=col-tile_cols_half, j2=0; j<col+tile_cols_half; j++, j2++)
      {
	if(i < 0 || j < 0 || i >= Q.rows() || j >= Q.cols())
	  result[i2][j2] = 0;
	else
	  result[i2][j2] = Q[i][j];
      }

  return result;
}

template<class T>
Appear<T> Train_Appear<T>::train(const _DMatrix<T> &part_locations, 
			   std::vector<std::string> &img_filenames,
			   int dilate, int tile_rows, int tile_cols)
{

  int part_count = part_locations.cols() / 2;


  assert(dilate == 5);
  
  const int dir_size = 16;
  T _bg_model[dir_size] =  { 0.4907, 0.0753, 0.0223, 0.0275, 0.0927,
				  0.0152,    0.0384,    0.0272, 
				  0.0229,    0.0272,    0.0070,    0.0228,    
				  0.0380, 0.0270, 	0.0283,    0.0376 };
				  
  _DMatrix<T> bg_model(1, dir_size, _bg_model);
  int image_count = part_locations.rows();
  
  _DMultiDMatrix<T> probs(4, part_count, dir_size, tile_rows, tile_cols);
  probs=0;

  for(int i=0; i<image_count; i++)
    {
      _DMatrix<char> Q = read_edge_map(img_filenames[i].c_str());

      
      for(int f=0; f<part_count; f++)
	{
	  _DMatrix<char> R=get_image_patch(Q, 
				    int(part_locations[i][f*2+1]), 
				    int(part_locations[i][f*2]), 
				    tile_rows, tile_cols);
	  
	  for(int row=0; row < tile_rows; row++)
	    for(int col = 0; col<tile_cols; col++)
	      {
		(probs.get(f, int(R[row][col]))[row][col])++;
	      }
	}
    }


  // probs = probs / image_count;

 for(int i=0; i<part_count; i++)
   for(int j=0; j<dir_size; j++)
     {
       _DMatrix<T> &matrix = probs.get(i,j);

       for(int k=0; k<tile_rows; k++)
	 for(int l=0; l<tile_cols; l++)
	   {
	     matrix[k][l] /= T(image_count);
	     if(matrix[k][l] == 0)
	       matrix[k][l] = 0.0001;
	   }
     }

 Appear<T> appear(probs, bg_model);

 return appear;
}

#define DECLARE(x) \
  template class Train_Appear<x>; 

DECLARE(double)
DECLARE(float)
