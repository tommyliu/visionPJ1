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
#include <vector>

template<class T>
class Train_Appear
{
 public:
  Train_Appear() {}

  Appear<T> train(const _DMatrix<T> &part_locations, 
	       std::vector<std::string> &img_filenames,
	       int dilate, int tile_rows, int tile_cols);

 protected:
  _DMatrix<char> read_edge_map(const char *fname);
  _DMatrix<char> get_image_patch(const _DMatrix<char> &Q, 
					  int row, int col, int tile_rows, 
					  int tile_cols);

};
