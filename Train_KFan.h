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
#include <KFan.h>

template<class T>
class Train_KFan
{
 public:
  Train_KFan() {};

  T glikelihood(const _DMatrix<T> &sigma, const _DMatrix<T> &mu, 
		     const _DMatrix<T> &data) const;

  KFan<T> train(const _DMatrix<T> &part_locations, const _DMatrix<int> &dist_parts);
  KFan<T> train(const _DMatrix<T> &part_locations, int dist_part_count);

 protected:



};
