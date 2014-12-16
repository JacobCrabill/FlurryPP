/*!
 * \file matrix.cpp
 * \brief Class for simplified matrix storage & manipulation
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

#include "global.hpp"


matrix::matrix()
{
  dim0 = 0;
  dim1 = 0;
}


matrix::matrix(int inDim0, int inDim1)
{
  data.resize(inDim0);
  for (auto& x: data) x.resize(inDim1);
}

void matrix::setup(int inDim0, int inDim1)
{
  data.resize(inDim0);
  for (auto& x: data) x.resize(inDim1);
}

T& matrix::operator[](int inDim0, int inDim1)
{
  if (inDim0<dim0 && inDim1<dim1) return data[inDim0][inDim1];
}


void matrix::initialize_to_zero()
{
  for (int idim=0; idim<dim0; idim++)
    for (int jdim=0; jdim<dim1; jdim++)
      data[idim][jdim] = 0;
}


template <typename T>
void matrix::timesMatrix(matrix<T> &A, matrix<T> &B)
{
  int i, j, k, p;
  p = A.dim1;

  if (A.dim0 != dim1) FatalError("Incompatible matrix sizes!");
  if (B.dim0 != dim0 || B.dim1 != A.dim1) B.setup(dim0, A.dim1);

  B.initialize_to_zero();

  for (i=0; i<dim0; i++) {
    for (j=0; j<dim1; j++) {
      for (k=0; k<p; k++) {
        B[i][k] += data[i][j]*A[j][k];
      }
    }
  }
}

template <typename T>
void matrix::timesVector(vector<T> &A, vector<T> &B)
{
  int i, j;

  if (A.size() != dim1) FatalError("Incompatible vector size");
  if (B.size() != dim1) B.resize(dim1);

  for (i=0; i<dim0; i++) {
    B[i] = 0;
    for (j=0; j<dim1; j++) {
      B[i] += data[i][j]*A[j];
    }
  }
}
