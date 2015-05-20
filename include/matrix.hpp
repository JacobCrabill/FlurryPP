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

  /* --- Member Functions --- */
  void setup(uint inDim0, uint inDim1);

  //! Adds the matrix a*A to current matrix (M += a*A)
  void addMatrix(matrix<T> &A, double a=1);

  //! Multiplies the matrix by the matrix A and stores the result in B (B = M*A)
  void timesMatrix(matrix<T> &A, matrix<T> &B);

  //! Multiplies the matrix by the matrix A and adds the result to B (B += M*A)
  void timesMatrixPlus(matrix<T> &A, matrix<T> &B);

  //! Multiplies the matrix by the vector A and stores the result in B (B = M*A)
  void timesVector(vector<T> &A, vector<T> &B);

  /*! Insert a row into the matrix at location rowNum [zero-indexed], with the default being at the end */
  void insertRow(vector<T> &vec, int rowNum = -1);

  void insertRow(T* vec, int rowNum, int length);

  void addCol(void);

  vector<T> getRow(uint row);

  matrix<T> getRows(vector<int> ind);

  vector<T> getCol(int col);

  void print(void);

  /* --- Data-Access Operators --- */

  T* operator[](int inDim0);

  T &operator()(int i, int j=0);

  T* getData();

  /* --- Search Operations --- */

  /*! Find all unique 'rows' in a matrix */
  void unique(matrix<T> &out, vector<int> &iRow);  

  /* --- Member Variables --- */
  uint dim0, dim1;  //! Dimensions of the matrix

protected:
  vector<T> data;
};
