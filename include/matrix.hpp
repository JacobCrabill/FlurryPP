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
  matrix<T>& operator=(const matrix<T>& inMatrix);

  //! Destructor
  ~matrix();

  /* --- Member Functions --- */
  void setup(int inDim0, int inDim1);

  //! Multiplies the oper by the matrix A and stores the result in B
  void timesMatrix(matrix<T> &A, matrix<T> &B);

  //! Multiplies the oper by the vector A and stores the result in B
  void timesVector(vector<T> &A, vector<T> &B);

  T& operator[](int inDim0, int inDim1);

  /* --- Member Variables --- */
  int dim0, dim1;  //! Dimensions of the matrix

private:
  vector<vector<T>> data;
};
