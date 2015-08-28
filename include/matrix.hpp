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
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#pragma once

#include <array>
#include <iomanip>   // for setw, setprecision
#include <iostream>
#include <vector>

#include "error.hpp"

struct point;

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

  //! Assignment
  Array<T,N> operator=(const Array<T,N>& inArray);

  /*! Get dim0 [number of rows] */
  uint getDim0(void) {return this->dims[0];}

  /*! Get dim1 [number of columns] */
  uint getDim1(void) {return this->dims[1];}

  /*! Get the size of the underlying data array (total number of Array elements) */
  uint getSize(void) {return data.size();}

  /* --- Member Functions --- */

  void setup(uint inDim0, uint inDim1=1, uint inDim2=1, uint inDim3=1);

  /* --- Data-Access Operators --- */

  /*! Returns a pointer to the first element of row inDim0 */
  T* operator[](int inDim0);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0, int k=0, int l=0);

  /*! Returns the .data() pointer of the underlying vector<T> data */
  T* getData();

  /* --- Member Variables --- */
  //uint dim0, dim1;  //! Dimensions of the Array

  array<uint,4> dims;  //! Dimensions of the Array

  vector<T> data;
};

template<typename T>
class Array2D : public Array<T,2>
{
public:

  Array2D();

  Array2D(uint inDim0, uint inDim1);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0);

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

  vector<T> getCol(int col);
};

template <typename T>
class matrix : public Array2D<T> {
public:

  matrix();

  matrix(uint inDim0, uint inDim1);

  void initializeToZero(void);

  void initializeToValue(T val);

  //! Adds the Array a*A to current Array (M += a*A)
  void addMatrix(matrix<T> &A, double a=1);

  //! Add another Array to the current Array
  matrix<double>& operator+=(matrix<double> &A);

  //! Subtract another Array to the current Array
  matrix<double>& operator-=(matrix<double> &A);

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

  /*! Invert the current Matrix using Gauss elimination with pivots */
  matrix<T> invertMatrix(void);

  //! Get the determinant of the matrix (recursive algorithm)
  double det();

  //! Get the adjoint of the matrix (transpose of the cofactor matrix)
  matrix<T> adjoint(void);

  /*! Reshape a 1D vector into a 2D Array */
  void vecToMatrixResize(vector<T> &A);

  /*! Prints the contents of the Array to the console */
  void print(void);

  /*! Check the matrix for NaN values */
  bool checkNan(void);

  /* --- Search Operations --- */

  /*! Find all unique 'rows' in a Array */
  void unique(matrix<T> &out, vector<int> &iRow);
};
