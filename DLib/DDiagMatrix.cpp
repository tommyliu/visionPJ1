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
#include "DDiagMatrix.h"
#include <string>

using namespace std;

DDiagMatrix::DDiagMatrix()
{
    _rows = _cols = 0;
}


DDiagMatrix::DDiagMatrix(int __rows, int __cols)
{
  _rows = __rows;
  _cols = __cols;
  
  data = DMatrix(1, _rows);
}


DDiagMatrix::DDiagMatrix(const DMatrix &source_data)
{
  assert(source_data.rows() == 1);
  
  data = source_data;
  _rows = source_data.cols();
  _cols = source_data.cols();
  
}

DDiagMatrix DDiagMatrix::operator+(const DDiagMatrix &other) const
{
  DMatrix result = data + other.data;
  
  return DDiagMatrix(result);
}

DDiagMatrix DDiagMatrix::operator+(double value) const
{
  DMatrix result = data + value;
  
  return DDiagMatrix(result);
}


DDiagMatrix DDiagMatrix::operator-(const DDiagMatrix &other) const
{
  DMatrix result = data - other.data;
  
  return DDiagMatrix(result);
}

DDiagMatrix DDiagMatrix::operator-(double value) const
{
  DMatrix result = data - value;
  
  return DDiagMatrix(result);
}

DDiagMatrix DDiagMatrix::operator*(const DDiagMatrix &other) const
{
  return DDiagMatrix(pointwise_multiply(data, other.data));
}

DDiagMatrix DDiagMatrix::operator*(double value) const
{
  return DDiagMatrix(data * value);
}


DDiagMatrix DDiagMatrix::operator/(double value) const
{
  return DDiagMatrix(data / value);
}


DDiagMatrix operator*(double value, const DDiagMatrix &other)
{
    return(other * value);
}


double DDiagMatrix::operator[](int entry) const
{
    return data[0][entry];
}

DDiagMatrix &DDiagMatrix::operator=(const DDiagMatrix &other)
{
    _rows = other.rows();
    _cols = other.cols();

    data = other.data;

    return *this;
}

DDiagMatrix::DDiagMatrix(const DDiagMatrix &other)
{
    *this = other;
}


DDiagMatrix operator+(double value, const DDiagMatrix &other)
{
    return(other + value);
}

DDiagMatrix operator-(double value, const DDiagMatrix &other)
{
  return DDiagMatrix(value - other.data);
}

DMatrix DDiagMatrix::operator*(const DMatrix &other) const
{
  assert(cols() == other.rows());

  DMatrix result(rows(), other.cols());

  for(int i=0; i<rows(); i++)
     for(int j=0; j<other.cols(); j++)
          result[i][j] = data[0][i] * other[i][j];

    return result;
}

void DDiagMatrix::set_entry(int row, int col, double value)
{
  assert(row == col);

  data[0][row] = value;
}

DMatrix operator*(const DMatrix &reg, const DDiagMatrix &diag) 
{
  assert(reg.cols() == diag.rows());

  DMatrix result(reg.rows(), diag.cols());

  for(int i=0; i<reg.rows(); i++)
    for(int j=0; j<diag.cols(); j++)
        result[i][j] = reg[i][j] * diag.data[0][j];

  return result;
}

DMatrix DDiagMatrix::operator+(const DMatrix &other) const
{
  assert(other.rows() == cols() && other.rows() == rows());

  DMatrix result(other.rows(), other.cols());
  result = other;

  for(int i=0; i<other.rows(); i++)
    result[i][i] += (*this)[i];

  return result;
}

DMatrix DDiagMatrix::operator-(const DMatrix &other) const 
{
  assert(other.rows() == other.cols() && other.rows() == rows());

  DMatrix result(other.rows(), other.cols());
  result = -other;

  for(int i=0; i<other.rows(); i++)
    result[i][i] += (*this)[i];

  return result;
}

ostream &operator<<(ostream &os, const DDiagMatrix &matrix)
{
    for(int i=0; i<matrix.rows(); i++)
    {
        for(int j=0; j<matrix.cols(); j++)
	  if(i != j)
	    os << "0 ";
	  else
	    os << matrix[i] << " ";
        os << endl;
    }

    return os;
}


DDiagMatrix::operator DMatrix() const
{
  DMatrix result(rows(), cols());

  for(int i=0; i<rows(); i++)
    for(int j=0; j<cols(); j++)
      if(i != j)
	result[i][j] = 0;
      else
	result[i][j] = data[0][i];

  return result;
}
