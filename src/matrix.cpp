/*
 * matrix.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Rajendra
 */

#ifndef MATRIX_CPP_
#define MATRIX_CPP_
#include "matrix.h"
#include <iostream>

using namespace std;
//Parameter Constructor
template<typename T>
Matrix<T>::Matrix(unsigned rows,unsigned cols,const T& initial){
	mat.resize(rows);
	for(unsigned i=0;i<mat.size();i++){
		mat[i].resize(cols,initial);
	}
	this->rows = rows;
	this->cols = cols;
}

//constructor from vector<vector>.
//As 2-D array passing in C++/C is pain, I rather choose the vector instead for initializing the matrix.
template<typename T>
Matrix<T>::Matrix(unsigned rows,unsigned cols,std::vector<std::vector<T> >& rhs){
	this->mat  = rhs;
	this->rows = rows;
	this->cols = cols;
}

//Constructor from vector.
template<typename T>
Matrix<T>::Matrix(unsigned rows,unsigned cols,std::vector<T>& rhs){
	std::vector<std::vector<T> > mat1;
	int vect_size = rhs.size();
	mat1.push_back(rhs);
	this->mat = mat1;
	this->rows = 1;
	this->cols = vect_size;
}

//Copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs){
	this->mat = rhs.mat;
	this->rows = rhs.get_rows();
	this->cols = rhs.get_cols();
}

//(Virtual) destructor
template<typename T>
Matrix<T>::~Matrix(){}

//Creation of Identity matrix..
//Looks like i don't need to define static again wierd c++ :)
template<typename T>
Matrix<T> Matrix<T>::eye(const unsigned& rhs){
	Matrix<T> result(rhs,rhs,0);
	for(unsigned i=0;i<rhs;i++){
		for(unsigned j=0;j<rhs;j++){
			if(i==j){
				result(i,j) = 1;
			}
		}
	}
	return result;
}


//Assignment Operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs){
	if(&rhs == this)
		return *this;

	unsigned new_rows = rhs.get_rows();
	unsigned new_cols = rhs.get_cols();

	mat.resize(new_rows);
	for(unsigned i=0;i<mat.size();i++){
		mat[i].resize(new_cols);
	}

	for(unsigned i=0;i<new_rows;i++){
		for(unsigned j=0;j<new_cols;j++){
			mat[i][j] = rhs(i,j);
		}
	}
	this->rows = new_rows;
	this->cols = new_cols;
	return *this;
}

//Addition of two matrices
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs){
	Matrix result(rows,cols,0.0);

	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result(i,j) = this->mat[i][j]+rhs(i,j);
		}
	}
	return result;
}

//Cumulative addition of this matrix and another
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs){
	unsigned rows = rhs.get_rows();
	unsigned cols = rhs.get_cols();
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			this->mat[i][j]+=rhs(i,j);
		}
	}
	return *this;
}

// Subtraction of this matrix and another
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) {
	unsigned rows = rhs.get_rows();
	unsigned cols = rhs.get_cols();
	Matrix result(rows,cols,0.0);
	for(unsigned i =0; i<rows ; i++) {
		for(unsigned j =0; j<cols ; j++) {
			result(i,j) = this->mat[i][j]- rhs(i,j);
		}
	}
	return result;
}

//Cumulative subtraction of this matrix and another
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs){
	unsigned rows = rhs.get_rows();
	unsigned cols = rhs.get_cols();
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			this->mat[i][j] -= rhs(i,j);
		}
	}
	return *this;
}

//Left multiplication of this matrix and other
//Shoudl work for non-square matrix too..
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs){

	unsigned rows_left = this->rows;
	unsigned cols_left = this->cols;
	unsigned rows_rhs = rhs.get_rows();
	unsigned cols_rhs = rhs.get_cols();

	//if(cols_left != rows_rhs){//throw error..}
	assert(cols_left == rows_rhs && "Matrix dimension must match..");

	Matrix result(rows_left,cols_rhs,0.0);
	for(unsigned i=0;i<rows_left;i++){
		for(unsigned j=0;j<cols_rhs;j++){
			double sum = 0;
			for(unsigned k=0;k<cols_left;k++){
				sum += this->mat[i][k] * rhs(k,j);
			}
			result(i,j) = sum;
		}
	}
	return result;
}

//Cumulative left multiplication of this matrix and another
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs){
	Matrix result = (*this)*rhs;
	(*this) = result;
	return *this;
}

//Calculate a transpose of this matrix
//Should be able to transpose no-square matrix too..
template<typename T>
Matrix<T> Matrix<T>::transpose(){
	Matrix result(cols,rows,0.0);
	for(unsigned i=0;i<cols;i++){
		for(unsigned j=0;j<rows;j++){
			result(i,j) = this->mat[j][i];
		}
	}
	return result;
}

//Calculate Minor matrix from parent matrix & calculate determinant and return
template<typename T>
Matrix<T> minor(Matrix<T>& rhs,int row,int col){

	//printMatrixx(rhs);
	//cout<<"i="<<row<<" j="<<col<<endl;
	int rows_ = rhs.get_rows();
	int cols_ = rhs.get_cols();

	Matrix<T> minor1(rows_-1,cols_-1,0);
	for(int i=0,ii=0;i<rows_;i++){
		if(i==row) continue;
		for(int j=0,jj=0;j<cols_;j++){
			if(j==col ) continue;
			minor1(ii,jj) = rhs(i,j);
			jj++;
		}
		ii++;
	}
	return minor1;
}

//Calculate pseudo inverse using moore-penrose method
//Looks like we need to do through svd instead..
/**
 * SVD of A is
 * A = U*(S 0;0 0)*V';
 * pinv(A) = V*(S^-1 0;0 0)*U';
 *
 */
template<typename T>
Matrix<T> Matrix<T>::pinverse(Matrix<T>& A){
	Matrix<T> pinv_A = A.transpose()*A;
	double det = determinant(pinv_A,pinv_A.get_rows());
	//printMatrixx(pinv_A);
	cout<<det<<endl;
	if( det != 0){
		pinv_A = inverse(pinv_A);
		pinv_A = pinv_A*A.transpose();
	}else{
		pinv_A = A*A.transpose();
		det = determinant(pinv_A,pinv_A.get_rows());
		printMatrixx(pinv_A);
		cout<<"dettt "<<det<<endl;
		pinv_A = inverse(pinv_A);
		pinv_A = A.transpose()*pinv_A;
	}

	return pinv_A;
}

//Calculate inverse of the matrix.
//You can follow this..
//http://timjones.tw/blog/archive/2014/10/20/the-matrix-inverted
template<typename T>
Matrix<T> Matrix<T>::inverse(Matrix<T>& rhs){
	//First need to check if it is non-singular or not i.e. determinant != 0
	int rows_ = rhs.get_rows();
	int cols_ = rhs.get_cols();

	assert(rows_== cols_ && "Matrix dimension must match i.e. square");
	T det = determinant(rhs,rows_);
	//cout<<"Matrix determinant = "<<det<<endl;
	assert(det != 0 && "Matrix should be non-singular i.e. determinant != 0 ");
	//If Matrix is just single element just return that..
	if(rows_ == 1 && cols_ ==1 )
		return rhs;
	//Find Minor matrix.
	Matrix<T> cofactor_mat(rows_,cols_,0);
	for(int i=0;i<rows_;i++){
		for(int j=0;j<cols_;j++){
			cofactor_mat(i,j) = pow(-1,(i+j))*determinant(minor(rhs,i,j),rows_-1);
		}
	}

	printMatrixx(cofactor_mat);
	Matrix<T> adjoint_mat = cofactor_mat.transpose();
	Matrix<T> inverse_mat = adjoint_mat/det;

	return inverse_mat;
}

//Calculate determinant of the matrix. Use recursion to calculate determinant..
template<typename T>
T Matrix<T>::determinant(Matrix<T> rhs,int n){
	//If dimension doesn't match throw error. Determinant is possible only for square matrix..i.e.e
	// rhs.get_row()== rhs.get_col()
	assert(rhs.get_rows()==rhs.get_cols() && "Matrix dimension must match i.e. matrix must be square..");
	T det = 0;
	int i,j,k,ii,jj;
	Matrix<double> sub_rhs(n-1,n-1,0);
	if(n==1){
		return rhs(0,0);
	}
	else if(n==2){
		return (rhs(0,0)*rhs(1,1)-rhs(0,1)*rhs(1,0));
	}else{
		for(k=0;k<n;k++){
			ii=0;
			for(i=1;i<n;i++){
				jj=0;
				for(j=0;j<n;j++){
					if(k==j){
						continue;
					}
					sub_rhs(ii,jj) = rhs(i,j);
					jj++;
				}
				ii++;
			}
			det = det + pow(-1,k)*rhs(0,k)*determinant(sub_rhs,n-1);
		}
	}
	return det;
}


//This function times is equivalent to .* in matlab. It is each element by element multiplication
//A=[1 1;2 2] B=[2 2;3 3] , times(A,B) = [1*2 1*2;2*3 2*3]
template<typename T>
Matrix<T> Matrix<T>::times(const Matrix<T>& rhs1, const Matrix<T>& rhs2){
	int row = rhs1.get_rows();
	int col = rhs1.get_cols();
	Matrix<T> result(row,col,0);
	for(unsigned i=0;i<row;i++){
		for(unsigned j=0;j<col;j++){
			result(i,j) = rhs1(i,j)*rhs2(i,j);
		}
	}
	return result;
}

//Append row to the given matrix. Where row is #total no. of row.
template<typename T>
Matrix<T>& Matrix<T>::appendRows(const unsigned& totalRow){
	int new_rows = this->rows + totalRow;
	Matrix<T> result(new_rows,cols,0);
	for(int i=0;i<new_rows;i++){
		for(int j=0;j<cols;j++){
			if(i<rows){
				result(i,j) = this->mat[i][j];
			}else{
				result(i,j) = 0;
			}
		}
	}
	(*this) = result;
	return *this;
}

//Append column to the given matrix. Where col is #total no. of column
template<typename T>
Matrix<T>& Matrix<T>::appendColumns(const unsigned& totalColumn){
	int new_cols = this->cols + totalColumn;
	Matrix<T> result(rows,new_cols,0);
	for(int i=0;i<rows;i++){
		for(int j=0;j<new_cols;j++){
			if(j<cols){
				result(i,j) = this->mat[i][j];
			}else{
				result(i,j) = 0;
			}
		}
	}
	(*this) = result;
	return *this;
}

//Will append matrix(rhs) on the existing(parent) matrix
template<typename T>
Matrix<T>& Matrix<T>::appendMatrix(const Matrix<T>& rhs){
	//Make sure the dimension of both matrix match. Need to add check/validation!! otherwise throw error..
	int new_rows = this->rows + rhs.get_rows();
	int new_cols = this->cols;
	int row_tmp = this->rows;

	Matrix<T> result(new_rows,new_cols,0);
	for(int i=0,k=0;i<new_rows;i++){
		if(i==row_tmp)k=0;
		for(int j=0;j<new_cols;j++){
			if(i<row_tmp){
				result(i,j) = this->mat[i][j];
			}
			else{
				result(i,j) = rhs(k,j);
			}
		}
		k++;
		cout<<endl;
	}
	(*this) = result;
	return *this;
}

//Will append column vector on parent columnNumber
template<typename T>
Matrix<T>& Matrix<T>::assignColumnVector(const Matrix<T>& columnVector, const unsigned& columnNumber){

	assert(this->rows == columnVector.get_rows() && "Parent column and vector column should be equal length");
	//assert(columnNumber <= total_col && "Column number should be 0 or greater and less than column size of the parent matrix");
	for(int i=0;i<this->rows;i++){
		this->mat[i][columnNumber] = columnVector(i,0);
	}
	return *this;
}
//Will append row vector on parent rowNumber
template<typename T>
Matrix<T>& Matrix<T>::assignRowVector(const Matrix<T>& rowVector, const unsigned& rowNumber){
	assert(this->cols == rowVector.get_cols() && "Parent row and vector row should be equal length");
	for(int j=0;j<this->cols;j++){
		this->mat[rowNumber][j] = rowVector(0,j);
	}
	return *this;
}

//Will delete column vector of m by n matrx and return m by n-1 matrix..
template<typename T>
Matrix<T>& Matrix<T>::deleteColumnVector(const unsigned& columnNumber){
	//Make sure rowNumber lies within the range..
	assert(columnNumber < cols && "columnNumber should lie within columns of parent matrix");
	Matrix<T> result(rows,cols-1,0);
	for(int i=0;i<rows;i++){
		for(int j=0,k=0;j<cols;j++){
			if(j == columnNumber)continue;
			result(i,k++) = this->mat[i][j];
		}
	}
	(*this) = result;
	return *this;
}

//Will delete row vector of m by n matrix and return m-1 by n matrix..
template<typename T>
Matrix<T>& Matrix<T>::deleteRowVector(const unsigned& rowNumber){
	//Make sure columnNumber lies within the range..
	assert(rowNumber < rows && "RowNumber to delete should lie in-between rows of parent matrix" );
	Matrix<T> result(rows-1,cols,0);
	for(int i=0,k=0;i<rows;i++,k++){
		if(i == rowNumber){
			k--;
			continue;
		}
		for(int j=0;j<cols;j++){
			result(k,j) = this->mat[i][j];
		}
	}
	(*this) = result;
	return *this;
}

//Return matrix from given row range. matrixRowRange(1,3) returns first 3 row of the matrix
template<typename T>
Matrix<T> Matrix<T>::matrixRowRange(unsigned fromRowId, unsigned toRowId){
	if(fromRowId==toRowId)
		toRowId = toRowId + 1;
	int new_row = toRowId - fromRowId;
	Matrix<T> result(new_row,cols,0);
	for(int i = fromRowId,k=0;i<toRowId;i++,k++){
		for(int j=0;j<cols;j++){
			result(k,j) = this->mat[i][j];
		}
	}
	return result;
}

//Return matrix from given column range. matrixColumnRange(2,3) returns 2 & 3 'rd row of the matrix
template<typename T>
Matrix<T> Matrix<T>::matrixColumnRange(unsigned fromColId, unsigned toColId){
	if(fromColId==toColId)toColId = toColId + 1;
	int new_col = toColId - fromColId;
	Matrix<T> result(rows,new_col,0);
	for(int i = 0;i<rows;i++){
		for(int j=fromColId,k=0;j<toColId;j++,k++){
			result(i,k) = this->mat[i][j];
		}
	}
	return result;
}


//Matrix/Scalar additon
template<typename T>
Matrix<T> Matrix<T>::operator+(const T& rhs){
	Matrix result(rows,cols,0.0);
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result(i,j) = this->mat[i][j]+rhs;
		}
	}
	return result;
}

//Matrix/scalar subtraction
template<typename T>
Matrix<T> Matrix<T>::operator-(const T& rhs){
	Matrix result(rows,cols,0.0);
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result(i,j) = this->mat[i][j] - rhs;
		}
	}
	return result;
}

//Matrix/scalar multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const T& rhs){
	Matrix result(rows,cols,0.0);
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result(i,j) = this->mat[i][j]*rhs;
		}
	}
	return result;
}

//Matrix/Scalar division
template<typename T>
Matrix<T> Matrix<T>::operator /(const T& rhs){
	Matrix result(rows,cols,0.0);
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result(i,j) = this->mat[i][j]/rhs;
		}
	}
	return result;
}


//Multiply a matrix with a vector
template<typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& rhs){
	std::vector<T> result(rhs.size(),0.0);
	for(unsigned i=0;i<rows;i++){
		for(unsigned j=0;j<cols;j++){
			result[i] = this->mat[i][j]*rhs[j];
		}
	}
	return result;
}

//Obtain a vector of the diagonal elements
template<typename T>
std::vector<T> Matrix<T>::diag_vec(){
	std::vector<T> result(rows,0.0);
	for(unsigned i=0;i<rows;i++){
		result[i]=this->mat[i][i];
	}
	return result;
}

//Get the Row vector from matrix specified by row number..
template<typename T>
Matrix<T> Matrix<T>::getRow(const unsigned& row){
	Matrix<T> matrix(1,cols,0.0);
	for(unsigned i=0;i<cols;i++){
		matrix(0,i) = this->mat[row][i];
	}
	return matrix;
}

//Get the column vector from matrix specified by column number..
template<typename T>
Matrix<T> Matrix<T>::getColumn(const unsigned& col){
	Matrix<T> matrix(rows,1,0.0);
	for(unsigned i=0;i<rows;i++){
		matrix(i,0) = this->mat[i][col];
	}
	return matrix;
}

//find the minimum value of the column given by columnId.
template<typename T>
T Matrix<T>::minColumn(const unsigned& colId){
	Matrix<T> matrix = getColumn(colId);
		T min = matrix(0,0);
		for(unsigned i = 0;i<rows;i++){
			if(matrix(i,0) < min){
				min = matrix(i,0);
			}
		}
	return min;
}

//find the maximum value of the column given by columnId.
template<typename T>
T Matrix<T>::maxColumn(const unsigned& colId){
	Matrix<T> matrix = getColumn(colId);
		T max = matrix(0,0);
		for(unsigned i = 0;i<rows;i++){
			if(matrix(i,0) > max){
				max = matrix(i,0);
			}
		}
		return max;
}

//find the maximum value of the row given by rowId.
template<typename T>
T Matrix<T>::maxRow(const unsigned& rowId){
	Matrix<T> matrix = getRow(rowId);
	T max = matrix(0,0);
	for(unsigned i = 0;i<cols;i++){
		if(matrix(0,i) > max){
			max = matrix(0,i);
		}
	}
	return max;
}

//find the minimum value of the row given by rowId..
template<typename T>
T Matrix<T>::minRow(const unsigned& rowId){
	Matrix<T> matrix = getRow(rowId);
	T min = matrix(0,0);
	for(unsigned i = 0;i<cols;i++){
		if(matrix(0,i) < min){
			min = matrix(0,i);
		}
	}
	return min;
}

//Access the individual elements
template<typename T>
T& Matrix<T>::operator()(const unsigned& row,const unsigned& col){
	return this->mat[row][col];
}

//Access the individual elements(const)
template<typename T>
const T& Matrix<T>::operator()(const unsigned& row,const unsigned& col) const{
	return this->mat[row][col];
}

//Get the number of rows of the matrix
template<typename T>
unsigned Matrix<T>::get_rows() const{
	return this->rows;
}

//Get the number of columns of the matrix
template<typename T>
unsigned Matrix<T>::get_cols() const{
	return this->cols;
}

//Prints the size of the matrix.
template<typename T>
void Matrix<T>::size() const{
	std::cout<<get_rows()<<" "<<get_cols()<<endl;
}

template<typename T>
void printMatrixx(Matrix<T> matt){
	for(int i=0;i<matt.get_rows();i++){
		for(int j=0;j<matt.get_cols();j++){
			cout<<matt(i,j)<<" ";
		}
		cout<<endl;
	}
}

#endif /* MATRIX_CPP_ */
