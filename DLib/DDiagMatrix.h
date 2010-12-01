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
#ifndef _DDIAGMATRIX_H_
#define _DDIAGMATRIX_H_

#include <iostream>
#include <DMatrix.h>
#include "DPoint.h"
#include "DRect.h"

class DDiagMatrix
{
    public:
        DDiagMatrix(int _rows, int _cols);
        DDiagMatrix();
        DDiagMatrix(const DMatrix &other);
	DDiagMatrix(const DDiagMatrix &other);

	operator DMatrix() const;

        ~DDiagMatrix() {}

        DDiagMatrix operator*(const DDiagMatrix &other) const;
        DDiagMatrix operator+(const DDiagMatrix &other) const;
        DDiagMatrix operator-(const DDiagMatrix &other) const;

        DDiagMatrix &operator=(const DDiagMatrix &other);

        DDiagMatrix operator+(double value) const;
        DDiagMatrix operator-(double value) const;
        DDiagMatrix operator*(double value) const;
        DDiagMatrix operator/(double value) const;

	DMatrix operator*(const DMatrix &other) const;
	DMatrix operator+(const DMatrix &other) const;
	DMatrix operator-(const DMatrix &other) const;
	friend DMatrix operator*(const DMatrix &reg, const DDiagMatrix &diag);

	double operator[](int index) const;
	void set_entry(int row, int col, double value);

        int rows() const { return _rows; }
        int cols() const { return _cols; }

        friend DDiagMatrix operator+(double value, const DDiagMatrix &other);
        friend DDiagMatrix operator-(double value, const DDiagMatrix &other);
        friend DDiagMatrix operator*(double value, const DDiagMatrix &other);

	friend std::ostream &operator<<(std::ostream &os, const DDiagMatrix &matrix);

    protected:
	DMatrix data;

        int _rows, _cols;
};


#endif
