/*
 * matrix.h
 *
 *  Created on: Feb 4, 2015
 *      Author: Rajendra
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <math.h>
#include <assert.h>
template<typename T> class Matrix{
private:
	std::vector<std::vector<T> > mat;
	unsigned rows;
	unsigned cols;

public:
	Matrix(unsigned rows,unsigned cols,const T& initial);
	Matrix(unsigned rows,unsigned cols,std::vector<std::vector<T> >& rhs);
	Matrix(unsigned rows,unsigned cols,std::vector<T>& rhs);
	Matrix(const Matrix<T>& rhs);
	virtual ~Matrix();

	//Creation of Identity matrix..
	static Matrix<T> eye(const unsigned& rhs);

	//Operator overloading for standard mathmatical matrix operations
	Matrix<T>& operator=(const Matrix<T>& rhs);

	//Matrix mathematical operations
	Matrix<T> operator+(const Matrix<T>& rhs);
	Matrix<T>& operator+=(const Matrix<T>& rhs);
	Matrix<T> operator-(const Matrix<T>& rhs);
	Matrix<T>& operator-=(const Matrix<T>& rhs);
	Matrix<T> operator*(const Matrix<T>& rhs);
	Matrix<T>& operator*=(const Matrix<T>& rhs);

	Matrix<T> transpose();
	static T determinant(Matrix<T> rhs,int n);
	static Matrix<T> inverse(Matrix<T>& rhs);
	static Matrix<T> pinverse(Matrix<T>& rhs);

	//This function times is equivalent to .* in matlab. It is each element by element multiplication
	static Matrix<T> times(const Matrix<T>& rhs1, const Matrix<T>& rhs2);

	//Will append rows or columns of the given matrix. Had to do due to matlab implementation :)
	//Will increase blank row or column accordingly..
	Matrix<T>& appendRows(const unsigned& totalRow);
	Matrix<T>& appendColumns(const unsigned& totalColumn);

	//Will append matrix on the existing matrix
	Matrix<T>& appendMatrix(const Matrix<T>& rhs);

	//Will append column vector on parent columnNumber
	Matrix<T>& assignColumnVector(const Matrix<T>& columnVector, const unsigned& columnNumber);
	Matrix<T>& assignRowVector(const Matrix<T>& columnVector, const unsigned& rowNumber);

	//Will delete column vector of m by n matrx and return m by n-1 matrix..
	Matrix<T>& deleteColumnVector(const unsigned& columnNumber);
	//Will delete row vector of m by n matrix and return m-1 by n matrix..
	Matrix<T>& deleteRowVector(const unsigned& rowNumber);

	//Return matrix from given row range. (1,3) returns first 3 row of the matrix
	Matrix<T> matrixRowRange(unsigned fromRowId, unsigned toRowId);
	Matrix<T> matrixColumnRange(unsigned fromColId, unsigned toColId);

	//Matrix/scalar operations
	Matrix<T> operator+(const T& rhs);
	Matrix<T> operator-(const T& rhs);
	Matrix<T> operator*(const T& rhs);
	Matrix<T> operator/(const T& rhs);

	//Matrix/vector operations
	std::vector<T> operator*(const std::vector<T>& rhs);
	std::vector<T> diag_vec();

	//Vector related operation on matrix.
	Matrix<T> getRow(const unsigned& row);
	Matrix<T> getColumn(const unsigned& col);
	T minRow(const unsigned& row);
	T maxRow(const unsigned& row);
	T minColumn(const unsigned& col);
	T maxColumn(const unsigned& col);

	//Access the individual elements
	T& operator()(const unsigned& row,const unsigned& col);
	const T& operator()(const unsigned& row,const unsigned& col) const;

	//Access the row and column sizes
	void size() const;
	unsigned get_rows() const;
	unsigned get_cols() const;

};

#include "matrix.cpp"

#endif /* MATRIX_H_ */
