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
#include <Train_KFan.h>

using namespace std;

template<class T>
T Train_KFan<T>::glikelihood(const _DMatrix<T> &sigma, const _DMatrix<T> &mu, 
			       const _DMatrix<T> &data) const
{
  assert(data.cols() == 2);

  T det_sigma = sigma[0][0] * sigma[1][1] - sigma[0][1] * sigma[1][0];
  T log_c = T( -(1/2.0) * data.cols() * data.rows() * log(2*M_PI) - 
    (1/2.0) * data.rows() * log((double)det_sigma));

  _DMatrix<T> sigma_inv(2,2);
  sigma_inv[0][0] = -sigma[1][1];
  sigma_inv[1][1] = -sigma[0][0];
  sigma_inv[0][1] = sigma[0][1];
  sigma_inv[1][0] = sigma[1][0];
  sigma_inv = sigma_inv / (-det_sigma);

  T e=0;
  for(int i=0; i<data.rows(); i++)
    {
      _DMatrix<T> _data = data.extract_row(i);
      e += ((_data - mu) * sigma_inv * (_data - mu).transpose() / 2.0)[0][0];
    }

  return( log_c - e);
}


template<class T>
KFan<T> Train_KFan<T>::train(const _DMatrix<T> &part_locations, int dist_part_count)
{
  switch(dist_part_count)
    {
    case 1:
      {
	KFan<T> best_kfan;
	_DMatrix<int> best_dist_parts;
	T best_l=-1e100;
	_DMatrix<int> dist_parts(1,1);
	
	for(int i=0; i<dist_part_count; i++)
	  {
	    dist_parts[0][0] = i;
	    KFan<T> kfan = train(part_locations, dist_parts);
	    
	    if(kfan.model_likelihood > best_l)
	      {
		best_l = kfan.model_likelihood;
		best_dist_parts = dist_parts;
		best_kfan = kfan;
	      }
	  }
	
	return best_kfan;
      }
    case 0:
      {
	KFan<T> result;
	result.parts_count = part_locations.cols()/2;
	result.K = 0;

	return result;
      }

    default:
      throw std::string("Train_KFan doesn't currently support that value of K");
    }

  return KFan<T>();
}


template<class T>
KFan<T> Train_KFan<T>::train(const _DMatrix<T> &part_locations, 
		       const _DMatrix<int> &dist_parts)
{
  int part_count = part_locations.cols() / 2;
  int img_count = part_locations.rows();
  int dist_part_count = dist_parts.cols();

  // first learn a gaussian for the distinguished parts
  _DMatrix<T> d_D(img_count, dist_part_count * 2);
  _DMatrix<T> d_D_diff(img_count, (dist_part_count-1)*2);

  for(int j=0; j<img_count; j++)
    {
      T *d_D_row = d_D[j];
      T *part_loc_row = part_locations[j];

      for(int i=0; i<dist_part_count; i++)
	{
	  d_D_row[i*2] = part_loc_row[dist_parts[0][i]*2];
	  d_D_row[i*2+1] = part_loc_row[dist_parts[0][i]*2+1];
	}

      if(dist_part_count > 1)
	{
	  T *d_D_diff_row = d_D_diff[j];
	  for(int i=1; i<dist_part_count; i++)
	    {
	      d_D_diff_row[i*2] = d_D_row[i*2] - d_D_row[0];
	      d_D_diff_row[i*2+1] = d_D_row[i*2+1] - d_D_row[1];
	    }
	}
    }

  T l2=0;
  _DMatrix<T> sigma_123, mu_123;
  if(dist_part_count > 1)
    {
      sigma_123 = d_D_diff.covariance();
      mu_123 = d_D_diff.means();
      
      l2 = glikelihood(sigma_123, mu_123, d_D_diff);
    }

  T l = l2;

  _DMultiDMatrix<T> sigma(3, part_count, 2, 2);
  _DMultiDMatrix<T> mu(3, part_count, 2, 1);
  _DMultiDMatrix<T> mu_x2(3, part_count, 2, (dist_part_count-1)*2);

  sigma=0;
  mu=0;
  mu_x2=0;


  // now learn gaussians for (p1, p2, ... , i) tuples for all other
  // parts i

  for(int i=0; i<part_count; i++)
    {
      bool d_part = false;
      for(int j=0; j<dist_part_count; j++)
	if(i == dist_parts[0][j])
	  d_part = true;

      if(d_part) continue;

      _DMatrix<T> this_D(img_count, dist_part_count * 2);
      for(int j=0; j<img_count; j++)
	{
	  //	  for(int k=1,k2=0; k<dist_part_count-1; k++,k2++)
	  //	    {
	  //	      this_D[j][k2*2] = d_D[j][k2*2] - d_D[j][0];
	  //	      this_D[j][k2*2+1] = d_D[j][k2*2+1] - d_D[j][1];
	      //	    }

	  //	  this_D[j][(dist_part_count-1)*2] = part_locations[j][i*2] - 
	  //	    d_D[j][0];
	  //	  this_D[j][(dist_part_count-1)*2+1] = part_locations[j][i*2+1] - 
	  //	    d_D[j][1];

	  //	  for(int k=0, k2=1; k<part_count-1; k++, k2++)
	  //	    {
	  int k2=i;
	  int k=0;
	  this_D[j][k*2] = part_locations[j][k2*2] - d_D[j][0];
	  this_D[j][k*2+1] = part_locations[j][k2*2+1] - d_D[j][1];
	      //	    }
	}

      _DMatrix<T> sigma_123i = this_D.covariance();

      _DMatrix<T> mu_123i = (this_D.means()).transpose();

      l=l+glikelihood(sigma_123i, mu_123i.transpose(), this_D) - l2;


      // now compute parameters for P(i|dist_parts), which is also normal
      int low_dist = 0, high_dist = (dist_part_count-1)*2-1;
      int low_i = high_dist+1, high_i = high_dist+2;

      _DMatrix<T> S_11 = sigma_123i.extract(DRect(low_i, low_i, high_i, high_i));

      if(dist_part_count != 1)
	{
	  _DMatrix<T> S_22 = sigma_123i.extract(DRect(low_dist, low_dist, high_dist, high_dist));
	  _DMatrix<T> S_12 = sigma_123i.extract(DRect(low_i, low_dist, high_i, high_dist));
	  _DMatrix<T> S_21 = S_12.transpose();

	  sigma.get(i) = S_11 - S_12 * S_22.inverse() * S_21;
	  mu_x2.get(i) = S_12 * S_22.inverse();
	  mu.get(i) = mu_123i.extract(DRect(low_i,0, high_i, mu_123i.cols()-1)) - 
	    S_12 * S_22.inverse() * 
	    mu_123i.extract(DRect(low_dist, 0, high_dist, mu_123i.cols()-1));
	}
      else
	{
	  sigma.get(i) = S_11;
	  mu.get(i) = mu_123i.extract(DRect(low_i,0, high_i, mu_123i.cols()-1));
	}


    }

 l = l / img_count;

 KFan<T> kfan(part_count, dist_part_count, dist_parts, sigma, mu, 
	   sigma_123, mu_123, mu_x2);
	   

 kfan.model_likelihood = l;

 return kfan;
}


#define DECLARE(x) \
  template class Train_KFan<x>; 

DECLARE(double)
DECLARE(float)




