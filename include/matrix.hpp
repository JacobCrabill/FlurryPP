/*!
 * \file matrix.hpp
 * \brief Header file for matrix class
 *
 * The matrix class is really just a handy wrapper for a vector<vector<T>>, with
 * some functions for setting up and multiplying
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

#include <vector>

using namespace std;

// Forward Declarations
template <typename T> class matrix;
template <typename T> class subMatrix;

template <typename T>
class matrix
{
public:
  /* --- Constructors, destructors --- */
  //! Default Constructor
  matrix();

  //! Secondary Constructor with Size Allocation
  matrix(int inDim0, int inDim1);

  //! Copy Constructor
  matrix(const matrix<T>& inMatrix);

  //! Assignment
  matrix<T> operator=(const matrix<T>& inMatrix);

  matrix<T> operator=(const subMatrix<T>& inSubMatrix);

  //! Destructor
  //~matrix();

  void initialize_to_zero();

  /*! Get dim0 [number of rows] */
  unsigned int getDim0(void) {return dim0;}

  /*! Get dim1 [number of columns] */
  unsigned int getDim1(void) {return dim1;}

  /* --- Member Functions --- */
  void setup(int inDim0, int inDim1);

  //! Multiplies the matrix by the matrix A and stores the result in B
  void timesMatrix(matrix<T> &A, matrix<T> &B);

  //! Multiplies the matrix by the vector A and stores the result in B
  void timesVector(vector<T> &A, vector<T> &B);

  /*! Insert a row into the matrix at location rowNum [zero-indexed], with the default being at the end */
  void insertRow(vector<T> &vec, int rowNum = -1);

  matrix<T> getRows(vector<int> ind);

  vector<T> getCol(int col);

  void print(void);

  /* --- Data-Access Operators --- */

  vector<T>& operator[](int inDim0);

  subMatrix<T> operator[](vector<int> &iRows);

  /* --- Search Operations --- */

  /*! Find all unique 'rows' in a matrix */
  void unique(matrix<T> &out, vector<int> &iRow);

  /* --- Member Variables --- */
  unsigned int dim0, dim1;  //! Dimensions of the matrix

  vector<vector<T> > getData();

//protected:

  vector<vector<T> > data;
};


template <typename T>
class subMatrix: public matrix<T>
{
public:
  subMatrix();

  subMatrix(matrix<T> *inMat, vector<int> iRows);

  subMatrix(matrix<T> *inMat, vector<int> &inRows, vector<int> &inCols);

  // have the matrix<t>::operator[] create a submatrix which, upon destruction,
  // puts its values back into the original matrix
  subMatrix<T> operator=(matrix<T>& inMatrix);

  subMatrix<T> operator=(subMatrix<T>& inSubMatrix);

  //~subMatrix();

private:
  matrix<T>* mat; // Pointer to parent matrix
  vector<int> rows, cols;
};
