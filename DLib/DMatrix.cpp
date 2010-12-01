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
#include "DMatrix.h"
#include <string>
#include <vector>

using namespace std;

template<class T>
void _DMatrix<T>::deallocate_storage()
{
  if(data)
    {
      delete[] data;
      delete[] data_area;
    }
}

template<class T>
void _DMatrix<T>::initialize_storage()
{
    if(data)
        deallocate_storage();

    data = new T *[_rows];
    data_area = new T[_rows * _cols];

    for(int i=0; i<_rows; i++)
        data[i] = &(data_area[i*_cols]);
}

template<class T>
_DMatrix<T>::_DMatrix()
{
    data = 0;
    data_area = 0;
    _rows = _cols = 0;
}

template<class T>
_DMatrix<T>::_DMatrix(int __rows, int __cols)
{
    _rows = __rows;
    _cols = __cols;

    data = 0;
    data_area = 0;

    initialize_storage();
}

template<class T>
_DMatrix<T>::_DMatrix(int __rows, int __cols, const T *array)
{
    _rows = __rows;
    _cols = __cols;

    data = 0;
    data_area = 0;

    initialize_storage();

    memcpy(data_area, array, _rows * _cols * sizeof(T));
}

template<class T>
bool same_size(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
    return(m1.rows() == m2.rows() && m1.cols() == m2.cols());
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator+(const _DMatrix<T> &other) const
{
    _DMatrix<T> result(_rows, _cols);

    if(!same_size(*this, other))
        throw string("Size mismatch in DMatrix operator +");

    T *cp1 = data_area, *cp2 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp2++, cp_out++)
            *cp_out = *cp1 + * cp2;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator+(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp_out++)
            *cp_out = *cp1 + value;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator-(const _DMatrix<T> &other) const
{
    _DMatrix<T> result(_rows, _cols);

    if(!same_size(*this, other))
        throw string("Size mismatch in DMatrix operator -");

    T *cp1 = data_area, *cp2 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; ++j, ++cp1, ++cp2, ++cp_out)
            *cp_out = *cp1 - * cp2;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator-(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; ++j, ++cp1, ++cp_out)
            *cp_out = *cp1 - value;

    return result;
}

template<class T>
_DMatrix<T> operator-(const _DMatrix<T> &m)
{
  _DMatrix<T> result(m.rows(), m.cols());

  T *cp_in = m[0], *cp_out = result[0];

  for(int i=0; i<m.rows(); ++i)
    for(int j=0; j<m.cols(); ++j, ++cp_in, ++cp_out)
      *cp_out = -*cp_in;

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator*(const _DMatrix<T> &other) const
{
    if(_cols != other._rows)
        throw string("Size mismatch in DMatrix operator *");

    _DMatrix<T> result(_rows, other._cols);

    for(int i=0; i<_rows; i++)
      {
	T *data_i = data[i];

        for(int j=0; j<other._cols; j++)
	  {
            T res=0;
	    
            for(int k=0; k<other._rows; k++)
	      res += data_i[k] * other.data[k][j];
	    
            result[i][j] = res;
	  }
      }

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator*(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp_out++)
            *cp_out = *cp1 * value;

    return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator/(T value) const
{
    _DMatrix<T> result(_rows, _cols);

    T *cp1 = data_area, *cp_out = result.data_area;
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++, cp1++, cp_out++)
            *cp_out = *cp1 / value;

    return result;
}

template<class T>
_DMatrix<T> operator*(T value, const _DMatrix<T> &other)
{
    return(other * value);
}

template<class T>
T *_DMatrix<T>::operator[](int row) const
{
    return data[row];
}

template<class T>
T _DMatrix<T>::operator[](DPoint &point) const
{
    return data[point.row()][point.col()];
}

template<class T>
_DMatrix<T>::~_DMatrix()
{
    deallocate_storage();
}

template<class T>
_DMatrix<T> &_DMatrix<T>::operator=(const _DMatrix<T> &other)
{
  assert(this != &other);
  
    _rows = other.rows();
    _cols = other.cols();

    //    data = 0;
    //    data_area = 0;

    initialize_storage();

    memcpy(data_area, other.data_area, _rows * _cols * sizeof(T));

    return *this;
}

template<class T>
_DMatrix<T> &_DMatrix<T>::operator=(T other)
{
    initialize_storage();

    T *cp = (*this)[0];
    for(int i=0; i<rows(); i++)
      for(int j=0; j<cols(); j++)
	*(cp++) = other;

    return *this;
}

template<class T>
_DMatrix<T>::_DMatrix(const _DMatrix<T> &other)
{
  assert(this != &other);
  
  data = 0;
  data_area = 0;

  *this = other;
}

template<class T>
_DMatrix<T> _DMatrix<T>::transpose() const
{
    _DMatrix<T> result(_cols, _rows);

    T *in_cp = data[0];
    for(int i=0; i<_rows; i++)
        for(int j=0; j<_cols; j++)
            result.data[j][i] = *(in_cp++);

    return result;
}

template<class T>
_DMatrix<T> operator+(T value, const _DMatrix<T> &other)
{
    return(other + value);
}

template<class T>
_DMatrix<T> operator-(T value, const _DMatrix<T> &other)
{
    _DMatrix<T> result(other.rows(), other.cols());

    T *cp1 = other.data_area, *cp_out = result.data_area;
    for(int i=0; i<other.rows(); i++)
        for(int j=0; j<other.cols(); j++, cp1++, cp_out++)
            *cp_out = -*cp1 + value;

    return result;
}

template<class T>
istream &operator>>(istream &is, _DMatrix<T> &matrix)
{
  is >> matrix._rows >> matrix._cols;

  matrix.data = 0;
  matrix.data_area = 0;

  matrix.initialize_storage();


  for(int i=0; i<matrix.rows(); i++)
    for(int j=0; j<matrix.cols(); j++)
      is >> matrix.data[i][j];
  
  return is;
}

template<class T>
ostream &operator<<(ostream &os, const _DMatrix<T> &matrix)
{
  os << matrix.rows() << " " << matrix.cols() << endl;

  for(int i=0; i<matrix.rows(); i++)
    {
      for(int j=0; j<matrix.cols(); j++)
	os << matrix.data[i][j] << " ";
      os << endl;
    }
  
  return os;
}

template<class T>
_DMatrix<T> _DMatrix<T>::LU_factor()
{
    _DMatrix<T> result((*this));

    for(int j=0; j<cols(); j++)
        for(int i=j+1; i<rows(); i++)
        {
            T alpha = result[j][j] / result[j][i];
            _DMatrix<T> this_row = extract_row(i);
            result.set_row(i, this_row-this_row*alpha);
        }
    
    return result; 
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_row(int row) const
{
    return extract(DRect(row, 0, row, cols()-1));
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract_col(int col) const
{
    return extract(DRect(0, col, rows()-1, col));
}

template<class T>
_DMatrix<T> _DMatrix<T>::extract(const DRect &rect) const
{
    assert(rect.top() >= 0 && rect.left() >= 0);
    assert(rect.right() < cols() && rect.bottom() < rows());

    _DMatrix<T> result(rect.height(), rect.width());

    T *cp = result[0];
    for(int i=rect.top(); i <= rect.bottom(); i++)
       for(int j=rect.left(); j <= rect.right(); j++, cp++)
           *cp = (*this)[i][j];

    return result;
}

template<class T>
void _DMatrix<T>::set_submatrix(const DPoint &pt, const _DMatrix<T> &in)
{
    assert(pt.row() + in.rows() <= rows());
    assert(pt.col() + in.cols() <= cols());

    for(int i=0; i<in.rows(); i++)
        for(int j=0; j<in.cols(); j++)
            (*this)[i+pt.row()][j+pt.col()]=in[i][j];
}

template<class T>
void _DMatrix<T>::set_row(int row, const _DMatrix<T> &in)
{
    set_submatrix(DPoint(row, 0), in);
}

template<class T>
void _DMatrix<T>::set_col(int col, const _DMatrix<T> &in)
{
    set_submatrix(DPoint(0, col), in);
}

template<class T>
_DMatrix<T> _DMatrix<T>::reshape(int new_rows, int new_cols)
{
    assert(new_rows * new_cols == rows() * cols());

    _DMatrix<T> result(new_rows, new_cols);

    T *in_cp = (*this)[0];
    T *out_cp = result[0];

    for(int i=0; i<new_rows*new_cols; i++, out_cp++, in_cp++)
       *out_cp = *in_cp;

    return result;
}

template<class T>
_DMatrix<T> sqrt(const _DMatrix<T> &m1)
{
  _DMatrix<T> result(m1.rows(), m1.cols());

  T *in_cp = m1[0], *out_cp = result[0];
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, in_cp++, out_cp++)
      *out_cp = sqrt(*in_cp);

  return result;
}

template<class T>
_DMatrix<T> pointwise_multiply(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++)
      result[i][j] = m1[i][j] * m2[i][j];

  return result;
}

template<class T>
_DMatrix<T> pointwise_divide(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++)
      result[i][j] = m1[i][j] / m2[i][j];

  return result;
}

template<class T>
_DMatrix<T> pointwise_min(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  T *cp1 = m1[0], *cp2 = m2[0], *cpr = result[0];
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, cp1++, cp2++, cpr++)
      if(*cp1 < *cp2)
	*cpr = *cp1;
      else
	*cpr = *cp2;

  return result;
}

template<class T>
_DMatrix<T> pointwise_max(const _DMatrix<T> &m1, const _DMatrix<T> &m2)
{
  assert(same_size(m1, m2));

  _DMatrix<T> result(m1.rows(), m2.cols());

  T *cp1 = m1[0], *cp2 = m2[0], *cpr = result[0];
  for(int i=0; i<m1.rows(); i++)
    for(int j=0; j<m1.cols(); j++, cp1++, cp2++, cpr++)
      if(*cp1 > *cp2)
	*cpr = *cp1;
      else
	*cpr = *cp2;

  return result;
}

template<class T>
_DMatrix<T> fabs(const _DMatrix<T> &m)
{
  _DMatrix<T> result(m.rows(), m.cols());

  for(int i=0; i<m.rows(); i++)
    for(int j=0; j<m.cols(); j++)
      result[i][j] = fabs(m[i][j]);

  return result;
}

// compute mean of all entries in matrix
template<class T>
T _DMatrix<T>::mean() const
{
  return sum() / (rows() * cols());
}


// compute sum of all entries in matrix
template<class T>
T _DMatrix<T>::sum() const
{
  T sum = 0;

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      sum += data[i][j];

  return sum;
}


// compute median of all entries in matrix
template<class T>
T _DMatrix<T>::median() const
{
  vector<T> vect;

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      vect.push_back((*this)[i][j]);

  sort(vect.begin(), vect.end());

  return vect[vect.size()/2];
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator==(T value)
{
  _DMatrix<T> result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      result[i][j] = (*this)[i][j] == value;

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator<(T value)
{
  _DMatrix<T> result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      result[i][j] = (*this)[i][j] < value;

  return result;
}

template<class T>
_DMatrix<T> _DMatrix<T>::operator>(T value)
{
  _DMatrix<T> result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      result[i][j] = (*this)[i][j] > value;

  return result;
}


// computes the mean of *each column* of the input matrix.
//
template<class T>
_DMatrix<T> _DMatrix<T>::means()
{
  _DMatrix<T> result(1, cols());

  for(int j=0; j<cols(); j++)
  {
    T mean=0;

    for(int i=0; i<rows(); i++)
      mean += (*this)[i][j];

    result[0][j] = mean / (rows());
  }

  return result;
}


// assumes each row is an observation, each column is
// a feature.
template<class T>
_DMatrix<T> _DMatrix<T>::covariance()
{
  _DMatrix<T> cov(cols(), cols());
  _DMatrix<T> mean_vector = means();
  T *mean_vec = mean_vector[0];

  for(int i=0; i<cols(); i++)
    for(int j=0; j<cols(); j++)
    {
       T var = 0;

       for(int k=0; k<rows(); k++)
	 var += ((*this)[k][i]-mean_vec[i]) * ((*this)[k][j]-mean_vec[i]);
 
       cov[i][j] = var / (rows()-1);
    }

  return cov;
}

template<class T>
void _DMatrix<T>::swap_rows(int row1, int row2)
{
  // ideally, we'd just have to swap the pointers, except that
  // some routines assume that the matrix is arranged
  // contiguously in row-major order. So we have to move
  // memory around also (unfortunately).

  T *temp = data[row1];
  data[row1] = data[row2];
  data[row2] = temp;

  temp = new T[cols()];
  memcpy(temp, data[row1], sizeof(T) * cols());
  memcpy(data[row1], data[row2], sizeof(T) * cols());
  memcpy(data[row2], temp, sizeof(T) * cols());

  delete[] temp;
}


#ifdef GSL_SUPPORT
template<class T>
gsl_matrix *DMatrix_to_gsl(_DMatrix<T> &dm)
{
  gsl_matrix *gm = gsl_matrix_alloc(dm.rows(), dm.cols());

  T *cp = dm[0];
  for(int i=0; i<dm.rows(); i++)
    for(int j=0; j<dm.cols(); j++, cp++)
      gsl_matrix_set(gm, i, j, *cp);

  return gm;
}

template<class T>
_DMatrix<T> gsl_to_DMatrix(gsl_matrix *gm)
{
  _DMatrix<T> dm(gm->size1, gm->size2);

  T *cp = dm[0];
  for(int i=0; i<dm.rows(); i++)
    for(int j=0; j<dm.cols(); j++, cp++)
      *cp = T(gsl_matrix_get(gm, i, j));

  return dm;
}

template<class T>
_DMatrix<T> gsl_to_DMatrix(gsl_vector *gm)
{
  _DMatrix<T> dm(1, gm->size);

  T *cp = dm[0];
  for(int i=0; i<dm.cols(); i++, cp++)
    *cp = T(gsl_vector_get(gm, i));

  return dm;
}

template<class T>
_DMatrix<T> _DMatrix<T>::inverse()
{
  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_permutation *p = gsl_permutation_alloc(rows());;
  int signum;

  gsl_linalg_LU_decomp(gm, p, &signum);

  gsl_matrix *inverse = gsl_matrix_alloc(rows(), cols());
  gsl_linalg_LU_invert(gm, p, inverse);

  _DMatrix<T> dm = gsl_to_DMatrix<T>(inverse);
  gsl_permutation_free(p);
  gsl_matrix_free(gm);
  gsl_matrix_free(inverse);

  return dm;
}

template<class T>
T _DMatrix<T>::determinant()
{
  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_permutation *p = gsl_permutation_alloc(rows());;
  int signum;

  gsl_linalg_LU_decomp(gm, p, &signum);
  gsl_matrix inverse;
  T det = T(gsl_linalg_LU_det(gm, signum));

  _DMatrix<T> dm = gsl_to_DMatrix<T>(gm);
  gsl_permutation_free(p);
  gsl_matrix_free(gm);

  return det;
}

template<class T>
std::pair<_DMatrix<T>, _DMatrix<T> > _DMatrix<T>::eigen()
{
  assert(rows() == cols());

  gsl_matrix *gm = DMatrix_to_gsl(*this);
  gsl_vector *eigval = gsl_vector_alloc(rows());
  gsl_matrix *eigvec = gsl_matrix_alloc(rows(), rows());
  gsl_eigen_symmv_workspace *worksp = gsl_eigen_symmv_alloc(rows());

  gsl_eigen_symmv(gm, eigval, eigvec, worksp);

  gsl_eigen_symmv_sort (eigval, eigvec, GSL_EIGEN_SORT_ABS_ASC);

  _DMatrix<T> d_eigvec = gsl_to_DMatrix<T>(eigvec);
  _DMatrix<T> d_eigval = gsl_to_DMatrix<T>(eigval);

  gsl_matrix_free(gm);
  gsl_matrix_free(eigvec);
  gsl_vector_free(eigval);
  gsl_eigen_symmv_free(worksp);

  return std::pair<_DMatrix<T>,_DMatrix<T> >(d_eigval, d_eigvec);

}




#endif


#define DECLARE(x) \
  template class _DMatrix<x>; \
  template _DMatrix<x> operator-(x value, const _DMatrix<x> &other); \
  template _DMatrix<x> pointwise_multiply(const _DMatrix<x> &m1, const _DMatrix<x> &m2); \
  template ostream &operator<<(ostream &os, const _DMatrix<x> &matrix); \
  template istream &operator>>(istream &is, _DMatrix<x> &matrix); \
  template _DMatrix<x> operator-(const _DMatrix<x> &m); \
  template _DMatrix<x> pointwise_min(const _DMatrix<x> &, const _DMatrix<x> &); \
  template _DMatrix<x> pointwise_max(const _DMatrix<x> &, const _DMatrix<x> &);



DECLARE(double)
DECLARE(short)
DECLARE(int)
DECLARE(float)
DECLARE(char)
DECLARE(unsigned char)
