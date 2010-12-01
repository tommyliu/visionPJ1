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
#ifndef __KFAN_H__
#define __KFAN_H__

#include <DImage.h>
#include <DMultiDMatrix.h>
#include <DBoxMinTransform.h>
#include <iostream>
#include <fstream>
#include <queue>

template<class T>
class KFan {

 public:
  KFan() {}
  KFan(int _parts, int _K, const _DMatrix<int> &_dist_parts, 
       const _DMultiDMatrix<T> &_cond_sigmas, 
       const _DMultiDMatrix<T> &_cond_means, 
       const _DMatrix<T> &_dist_sigmas, 
       const _DMatrix<T> &_dist_means,
       const _DMultiDMatrix<T> &_cond_means_x2)
    {
      parts_count = _parts;
      K = _K;
      dist_parts = _dist_parts;
      cond_sigmas = _cond_sigmas, cond_means = _cond_means;
      dist_sigmas = _dist_sigmas, dist_means = _dist_means;
      cond_mus_x2 = _cond_means_x2;
    }

  KFan(char *fname)
    {
      load_params(fname);
    }

  virtual std::pair<T, _DMatrix<T> > do_matching(_DMultiDMatrix<T> &_cost_maps) {}

  template<class T2>
  friend std::ostream &operator<<(std::ostream & ofs, KFan<T2> &kfan);
  
  void load_params(char *fname);

  // protected:
  int parts_count, K;
  _DMatrix<int> dist_parts;

  _DMultiDMatrix<T> cond_sigmas, cond_means, cond_mus_x2;
  _DMatrix<T> dist_sigmas, dist_means;
  T model_likelihood;
};

template<class T>
class K_1Fan : public KFan<T>
{
 public:
  K_1Fan(char *fname) : KFan<T>(fname) {dist_part=int(KFan<T>::dist_parts[0][0]);}
    
    virtual std::pair<T, _DMatrix<T> > do_matching(_DMultiDMatrix<T> &_cost_maps);

    int dist_part;
};

template<class T>
class K_0Fan : public KFan<T>
{
 public:
  K_0Fan(char *fname) : KFan<T>(fname) {}
    
    virtual std::pair<T, _DMatrix<T> > do_matching(_DMultiDMatrix<T> &_cost_maps);

    int dist_part;
};


#endif
