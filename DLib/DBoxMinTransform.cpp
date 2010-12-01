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
#include <DBoxMinTransform.h>
#include <math.h>

using namespace std;

template<class T>
void DBoxMinTransform<T>::del_prime_d(T *delta, int cols, int inc, 
				      int log_d, T **del_primes_2n)
{
  T *del_primes_2n_0 = del_primes_2n[0];
  if(inc == 1)
    memcpy(del_primes_2n_0, delta, sizeof(T) * cols);
  else
    for(int i=0, cp=0; i<cols; ++i, cp+=inc)
      del_primes_2n_0[i] = delta[cp];
 
  
  int n_minus1_size = 1;
  for(int n=1; n<=log_d; ++n, n_minus1_size*=2)
    {
      T *del_primes_2n_nminus1 = del_primes_2n[n-1];
      T *del_primes_2n_n = del_primes_2n[n];

      for(int i=0; i<cols-n_minus1_size; ++i)
	del_primes_2n_n[i] = MIN(del_primes_2n_nminus1[i],
				 del_primes_2n_nminus1[i + n_minus1_size]);

      for(int i=MAX(cols-n_minus1_size,0); i < cols; ++i)
	del_primes_2n_n[i] = (del_primes_2n_nminus1[i]);


    }

}

template<class T>
void DBoxMinTransform<T>::del_prime_d2(T *delta, int cols, int d, int inc,
				       int log_d, T **del_primes_2n)
{
  T result[cols];
  int n_size = 1;
  for(int i=0; i<cols; i++)
    result[i] = T(1e100);
  
  int sofar = 0;
  for(int n=0; n<=log_d; ++n, n_size *= 2)
    {
      if((d & n_size))
	{
	  T *del_primes_2n_n = del_primes_2n[n];

	  //	  assert(cols-sofar >= 0);
	  
	  for(int i=0; i<cols-sofar; ++i)
	    result[i] = MIN(result[i], del_primes_2n_n[i+sofar]);
	  
	  sofar += n_size;
	}
    }

  if(inc == 1)
    memcpy(delta, result, sizeof(T) * cols);
  else
    for(int i=0, cp=0; i<cols; ++i, cp+=inc)
      delta[cp] = result[i];
}


template<class T>
_DMatrix<T> DBoxMinTransform<T>::do_transform(int d_row, int d_col)
{
  int log_d = int(ceil(log2(double(d_col+1))));

  assert(log_d <= max_log_d);
  for(int i=0; i<Delta1.rows(); i++)
    {
      //      del_prime_d(Delta1[i], Delta.cols(), d_col+1, 1, log_d);
      del_prime_d2(Delta1[i], Delta1.cols(),  d_col+1, 1, log_d, del_primes_2n[i]);
    }

  log_d = int(ceil(log2(double(d_row+1))));
  assert(log_d <= max_log_d);
  T* tmp[log_d+1];
  for(int i=0; i<log_d+1; i++)
    tmp[i] = new T[Delta1.rows()];

  for(int i=0; i<Delta1.cols(); i++)
    {
      del_prime_d(Delta1[0]+i, Delta1.rows(), Delta1.cols(), log_d, tmp);
      del_prime_d2(Delta1[0]+i, Delta1.rows(), d_row+1, Delta1.cols(), log_d, tmp);
    }

  for(int i=0; i<log_d+1; i++)
    delete[] tmp[i];


  return Delta1;
}




template<class T>
DBoxMinTransform<T>::DBoxMinTransform(const _DMatrix<T> &img, int max_rows, int max_cols)
{
  Delta1 = img;
  max_log_d = MAX(int(ceil(log2(double(max_rows+1)))),
		      int(ceil(log2(double(max_cols+1)))));

  del_primes_2n = new T**[img.rows()];
  for(int i=0; i<img.rows(); i++)
    {
      del_primes_2n[i] = new T*[max_log_d+1];

      for(int j=0; j<max_log_d+1; j++)
	del_primes_2n[i][j] = new T[img.cols()];
    }

  //  std::cout << "------------------- " << del_primes_2n[84][2][0] << std::endl;

  for(int i=0; i<Delta1.rows(); i++)
    del_prime_d(Delta1[i], Delta1.cols(), 1, max_log_d, del_primes_2n[i]);
}



template<class T>
DBoxMinTransform<T>::~DBoxMinTransform()
{
  if(del_primes_2n)
    {
      for(int j=0; j<Delta1.rows(); j++)
	{
	  for(int i=0; i<max_log_d+1; i++)
	    delete[] del_primes_2n[j][i];

	  delete[] del_primes_2n[j];
	}
      
      delete[] del_primes_2n;
    }
}




template class DBoxMinTransform<double>;
template class DBoxMinTransform<float>;
template class DBoxMinTransform<int>;
template class DBoxMinTransform<short>;
template class DBoxMinTransform<char>;
