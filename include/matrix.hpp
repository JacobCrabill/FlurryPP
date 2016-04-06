/*!
 * \file Array.hpp
 * \brief Header file for Array class
 *
 * The Array class is really just a handy wrapper for a vector<vector<T>>, with
 * some functions for setting up and multiplying
 *
 * Yes, I'm well aware that I really should just be using boost::multi_array
 * instead for loads of reasons, but this was more fun.  Don't judge me.
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <array>
#include <iomanip>   // for setw, setprecision
#include <iostream>
#include <vector>

#include "error.hpp"

#include "global.hpp"

#define INSERT_AT_END -1

using namespace std;

typedef unsigned int uint;

template <typename T, uint N>
class Array
{
public:
  /* --- Constructors, destructors --- */
  //! Default Constructor
  Array();

  //! Secondary Constructor with Size Allocation
  Array(uint inDim0, uint inDim1=1, uint inDim2=1, uint inDim3=1);

  //! Copy Constructor
  Array(const Array<T,N>& inArray);

  void initializeToValue(const T &_val);

  //! Assignment
  Array<T,N> operator=(const Array<T,N>& inArray);

  /*! Get dim0 [number of rows] */
  uint getDim0(void) {return this->dims[0];}
  uint getDim0(void) const {return this->dims[0];}

  /*! Get dim1 [number of columns] */
  uint getDim1(void) {return this->dims[1];}
  uint getDim1(void) const {return this->dims[1];}

  /*! Get the size of the underlying data array (total number of Array elements) */
  uint getSize(void) {return data.size();}

  /* --- Member Functions --- */

  void setup(uint inDim0, uint inDim1=1, uint inDim2=1, uint inDim3=1);

  /* --- Data-Access Operators --- */

  /*! Returns a pointer to the first element of row inDim0 */
  T* operator[](int inDim0);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0, int k=0, int l=0);

  T operator()(int i, int j=0, int k=0, int l=0) const;

  /*! Returns the .data() pointer of the underlying vector<T> data */
  T* getData();

  /* --- Member Variables --- */
  //uint dim0, dim1;  //! Dimensions of the Array

  uint dims[4];  //! Dimensions of the Array

  vector<T> data;

  void add_dim_0(uint ind, const T &val);
  void add_dim_1(uint ind, const T &val);
  void add_dim_2(uint ind, const T &val);
  void add_dim_3(uint ind, const T &val);

  void remove_dim_0(uint ind);
  void remove_dim_1(uint ind);
  void remove_dim_2(uint ind);
  void remove_dim_3(uint ind);

};

template<typename T>
class Array2D : public Array<T,2>
{
public:

  Array2D();

  Array2D(uint inDim0, uint inDim1);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0);

  T operator()(int i, int j=0) const;

  /*! Append the input Array as rows at the end of the current Array */
  void appendRows(Array2D<T> &mat);

  /*! Insert a row into the Array at location rowNum [zero-indexed], with the default being at the end */
  void insertRow(const vector<T> &vec, int rowNum = -1);

  void insertRow(T* vec, int rowNum, int length);

  void insertRowUnsized(const vector<T> &vec);

  void insertRowUnsized(T* vec, uint length);

  //! Insert a column at the end of the Array
  void addCol(void);

  //! Insert a column at the end of the Array
  void addCols(int nCols);

  //! Remove columns from the end of the Array
  void removeCols(int nCols = 1);

  vector<T> getRow(uint row);

  //! Return a sub-Array view of the Array using the rows in ind
  Array2D<T> getRows(vector<int> ind);

  Array2D<T> slice(array<int, 2> rows, array<int, 2> cols);

  vector<T> getCol(int col);

  Array2D<T> transpose(void);
};

template <typename T>
class matrix : public Array2D<T> {
public:

  matrix();

  matrix(uint inDim0, uint inDim1);

  matrix(const Array2D<T> &inMatrix);

  void initializeToZero(void);

  //void initializeToValue(T val);

  //! Adds the Array a*A to current Array (M += a*A)
  void addMatrix(matrix<T> &A, double a=1);

  //! Add another Array to the current Array
  matrix<double>& operator+=(matrix<double> &A);

  //! Subtract another Array to the current Array
  matrix<double>& operator-=(matrix<double> &A);

  //! Multiply by scalar
  matrix<double>& operator*=(double a);

  //! Divide by scalar
  matrix<double>& operator/=(double a);

  //! Add two matrices
  friend matrix<T> operator+(matrix<T> A, const matrix<T> &B) { return A += B; }

  //! Subtract two matrices
  friend matrix<T> operator-(matrix<T> A, const matrix<T> &B) { return A -= B; }

  //! Multiplies the Array by the Array A and stores the result in B (B = M*A)
  void timesMatrix(matrix<T> &A, matrix<T> &B);

  //! Multiply two matrices
  //friend matrix<T> operator*(matrix<T> &A, matrix<T> &B);

  //! Multiplies the Array by the Array A and adds the result to B (B += M*A)
  void timesMatrixPlus(matrix<T> &A, matrix<T> &B);

  //! Multiplies the Array by the vector A and stores the result in B (B = M*A)
  void timesVector(vector<T> &A, vector<T> &B);

  /*! Invert the current matrix using Gauss elimination with pivots */
  matrix<T> invertMatrix(void);

  /*! Solve Ax=b with this matrix, given the RHS vector */
  vector<T> solve(vector<T> RHS);

  //! Get the determinant of the matrix (recursive algorithm)
  double det();

  //! Get the adjoint of the matrix (transpose of the cofactor matrix)
  matrix<T> adjoint(void);

  /*! Reshape a 1D vector into a 2D Array */
  void vecToMatrixResize(vector<T> &A);

  /*! Prints the contents of the Array to the console */
  void print(int prec=8) const;

  /*! Check the matrix for NaN values */
  bool checkNan(void);

  /*! Calculate Frobenius norm */
  T frobNorm(void);

  /* --- Search Operations --- */

  /*! Find all unique 'rows' in a Array */
  void unique(matrix<T> &out, vector<int> &iRow);
};
