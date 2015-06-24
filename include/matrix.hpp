/*!
 * \file matrix.hpp
 * \brief Header file for matrix class
 *
 * The matrix class is really just a handy wrapper for a vector<vector<T>>, with
 * some functions for setting up and multiplying
 *
 * Yes, I'm well aware that I really should just be using boost::multi_array
 * instead for loads of reasons, but this was more fun.  Don't judge me.
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#pragma once

#include <iomanip>   // for setw, setprecision
#include <iostream>
#include <vector>

#include "error.hpp"


#define INSERT_AT_END -1

using namespace std;

typedef unsigned int uint;

template <typename T>
class matrix
{
public:
  /* --- Constructors, destructors --- */
  //! Default Constructor
  matrix();

  //! Secondary Constructor with Size Allocation
  matrix(uint inDim0, uint inDim1);

  //! Copy Constructor
  matrix(const matrix<T>& inMatrix);

  //! Assignment
  matrix<T> operator=(const matrix<T>& inMatrix);

  void initializeToZero(void);

  void initializeToValue(T val);

  /*! Get dim0 [number of rows] */
  uint getDim0(void) {return dim0;}

  /*! Get dim1 [number of columns] */
  uint getDim1(void) {return dim1;}

  /*! Get the size of the underlying data array (total number of matrix elements) */
  uint getSize(void) {return data.size();}

  /* --- Member Functions --- */

  void setup(uint inDim0, uint inDim1);

  //! Adds the matrix a*A to current matrix (M += a*A)
  void addMatrix(matrix<T> &A, double a=1);

  //! Add another matrix to the current matrix
  matrix<double>& operator+=(matrix<double> &A);

  //! Subtract another matrix to the current matrix
  matrix<double>& operator-=(matrix<double> &A);

  //! Add two matrices
  friend matrix<T> operator+(matrix<T> A, const matrix<T> &B) { return A += B; }

  //! Subtract two matrices
  friend matrix<T> operator-(matrix<T> A, const matrix<T> &B) { return A -= B; }

  //! Multiplies the matrix by the matrix A and stores the result in B (B = M*A)
  void timesMatrix(matrix<T> &A, matrix<T> &B);

  //! Multiply two matrices
  //friend matrix<T> operator*(matrix<T> &A, matrix<T> &B);

  //! Multiplies the matrix by the matrix A and adds the result to B (B += M*A)
  void timesMatrixPlus(matrix<T> &A, matrix<T> &B);

  //! Multiplies the matrix by the vector A and stores the result in B (B = M*A)
  void timesVector(vector<T> &A, vector<T> &B);

  /*! Insert a row into the matrix at location rowNum [zero-indexed], with the default being at the end */
  void insertRow(const vector<T> &vec, int rowNum = -1);

  void insertRowUnsized(const vector<T> &vec);

  void insertRow(T* vec, int rowNum, int length);

  //! Insert a column at the end of the matrix
  void addCol(void);

  //! Insert a column at the end of the matrix
  void addCols(int nCols);

  //! Remove columns from the end of the matrix
  void removeCols(int nCols);

  vector<T> getRow(uint row);

  //! Return a sub-matrix view of the matrix using the rows in ind
  matrix<T> getRows(vector<int> ind);

  vector<T> getCol(int col);

  /*! Prints the contents of the matrix to the console */
  void print(void);

  /*! Invert the current matrix using Gauss elimination with pivots */
  matrix<T> invertMatrix(void);

  /*! Reshape a 1D vector into a 2D matrix */
  void vecToMatrixResize(vector<T> &A);

  /* --- Data-Access Operators --- */

  /*! Returns a pointer to the first element of row inDim0 */
  T* operator[](int inDim0);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0);

  /*! Returns the .data() pointer of the underlying vector<T> data */
  T* getData();

  /* --- Search Operations --- */

  /*! Find all unique 'rows' in a matrix */
  void unique(matrix<T> &out, vector<int> &iRow);  

  /* --- Member Variables --- */
  uint dim0, dim1;  //! Dimensions of the matrix

protected:
  vector<T> data;
};
