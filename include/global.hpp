/*!
 * \file global.hpp
 * \brief Header file for global constants, objects, and variables
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

#include <limits.h>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstddef>    // std::size_t
#include <cstdlib>
#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>

#ifdef _OMP
#include <omp.h>
#ifdef _MKL_BLAS
#include "mkl_types.h"
#include "mkl_cblas.h"
#else
#include "cblas.h"
#endif
#endif

#include "error.hpp"

template<typename T> class matrix;

#include "matrix.hpp"

// Forward declarations of basic Flurry classes
class geo;
class ele;
class solver;

using namespace std;

typedef unsigned int uint;

/* --- Misc. Common Constants / Globally-Useful Variables --- */

static double pi = 4.0*atan(1);

extern map<string,int> bcStr2Num;
//extern map<int,string> bcNum2Str;

/*! enumeration for element type */
enum ETYPE {
  TRI     = 0,
  QUAD    = 1,
  TET     = 2,
  PRISM   = 3,
  HEX     = 4,
  PYRAMID = 5
};

/*! Enumeration for original, mesh-file-defined face type */
enum FACE_TYPE {
  HOLE_FACE = -1,
  INTERNAL  = 0,
  BOUNDARY  = 1,
  MPI_FACE  = 2,
  OVER_FACE = 3
};

/*! Enumeration for original, mesh-file-defined node type */
enum NODE_TYPE {
  NORMAL_NODE = 0,
  OVERSET_NODE = 1,
  BOUNDARY_NODE = 2
};

/*! Enumeration for mesh (either create cartesian mesh or read from file) */
enum meshType {
  READ_MESH   = 0,
  CREATE_MESH = 1,
  OVERSET_MESH = 2
};

enum EQUATION {
  ADVECTION_DIFFUSION = 0,
  NAVIER_STOKES       = 1
};

/*! Enumeration for all available boundary conditions */
enum BC_TYPE {
  NONE = -1,
  PERIODIC = 0,
  CHAR_INOUT = 1,
  SUP_IN = 2,
  SUP_OUT = 3,
  SUB_IN = 4,
  SUB_OUT = 5,
  SUB_IN_CHAR = 6,
  SUB_OUT_CHAR = 7,
  SLIP_WALL = 8,
  ISOTHERMAL_NOSLIP = 9,
  ADIABATIC_NOSLIP = 10,
  OVERSET = 11,
  SYMMETRY = 12,
  ISOTHERMAL_NOSLIP_MOVING = 13
};

/*! Enumeration for VCJH scheme (correction function) to use */
enum VCJH_SCHEME {
  DG = 0,
  SD = 1,
  HU = 2,
  CPLUS = 3
};

/*! For convinience with geometry, a simple struct to hold an x,y,z coordinate */
struct point
{
  double x, y, z;

  point() {
    x = 0;
    y = 0;
    z = 0;
  }

  point (double _x, double _y, double _z) {
    x = _x;
    y = _y;
    z = _z;
  }

  point(double* pt, int nDims=3) {
    x = pt[0];
    y = pt[1];
    if (nDims==3)
      z = pt[2];
    else
      z = 0;
  }

  void zero() {
    x = 0;
    y = 0;
    z = 0;
  }

  double& operator[](int ind) {
    switch(ind) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        cout << "ind = " << ind << ": " << flush;
        FatalError("Invalid index for point struct.");
    }
  }

  double operator[](int ind) const {
    switch(ind) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        cout << "ind = " << ind << ": " << flush;
        FatalError("Invalid index for point struct.");
    }
  }

  point operator=(double* a) {
    struct point pt;
    pt.x = a[0];
    pt.y = a[1];
    pt.z = a[2];
    return pt;
  }

  point operator-(point b) {
    struct point c;
    c.x = x - b.x;
    c.y = y - b.y;
    c.z = z - b.z;
    return c;
  }

  point operator+(point b) {
    struct point c;
    c.x = x + b.x;
    c.y = y + b.y;
    c.z = z + b.z;
    return c;
  }

  point operator/(point b) {
    struct point c;
    c.x = x / b.x;
    c.y = y / b.y;
    c.z = z / b.z;
    return c;
  }

  point& operator+=(point b) {
    x += b.x;
    y += b.y;
    z += b.z;
    return *this;
  }

  point& operator-=(point b) {
    x -= b.x;
    y -= b.y;
    z -= b.z;
    return *this;
  }

  point& operator+=(double* b) {
    x += b[0];
    y += b[1];
    z += b[2];
    return *this;
  }

  point& operator-=(double* b) {
    x -= b[0];
    y -= b[1];
    z -= b[2];
    return *this;
  }

  point& operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }

  point& operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  double operator*(point b) {
    return x*b.x + y*b.y + z*b.z;
  }

  void abs(void) {
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);
  }

  double norm(void) {
    return std::sqrt(x*x+y*y+z*z);
  }

  point cross(point b) {
    point v;
    v.z = x*b.y - y*b.x;
    v.y = x*b.z - z*b.x;
    v.x = y*b.z - z*b.y;
    return v;
  }

};

point operator/(point a, double b);
point operator*(point a, double b);

bool operator<(const point &a, const point &b); // Just a sort of 'dummy' function for sorting purposes

std::ostream& operator<<(std::ostream &os, const point &pt);

std::ostream& operator<<(std::ostream &os, const matrix<double> &mat);
std::ostream& operator<<(std::ostream &os, const matrix<int> &mat);

std::ostream& operator<<(std::ostream &os, const vector<int> &vec);
std::ostream& operator<<(std::ostream &os, const vector<double> &vec);

double getDist(point a, point b);

//! For clearer notation when a vector is implied, rather than a point
typedef struct point Vec3;

matrix<double> createMatrix(vector<point> &pts);

int factorial(int n);

bool checkNaN(vector<double> &vec);

bool checkNaN(double* vec, int size);

/*! Get polynomial-order-based CFL limit.  Borrowed from Josh's zefr code. */
double getCFLLimit(int order);

/*! Find indices of all values in vec equal to val */
template<typename T>
vector<int> findEq(const vector<T> &vec, T val)
{
  vector<int> out;

  for (uint i=0; i<vec.size(); i++) {
    if (vec[i]==val) out.push_back(i);
  }

  return out;
}

/*! Find index of first occurance of val in vec */
template<typename T>
int findFirst(vector<T> &vec, T val)
{
  if (vec.size()==0) return -1;

  for (int i=0; i<(int)vec.size(); i++) {
    if (vec[i]==val) return i;
  }

  // If not found...
  return -1;
}

template<typename T>
int findFirst(T* vec, T val, uint length)
{
  if (length==0) return -1;

  for (int i=0; i<(int)length; i++) {
    if (vec[i]==val) return i;
  }

  // If not found...
  return -1;
}

/*! Assign a value to vector vec at indices indicated in ind */
template<typename T>
void vecAssign(vector<T> &vec, vector<int> &ind, T val)
{
  for (auto& i:ind) vec[i] = val;
}

// Good for numeric types - meant for ints or uints though
template<typename T>
T getMax(vector<T> &vec)
{
  T max = 0;
  for (auto& i:vec) {
    if (i>max) max = i;
  }

  return max;
}

// Good for numeric types
template<typename T>
T getMin(vector<T> &vec)
{
  T min = 1e15;
  for (auto& i:vec) {
    if (i<min) min = i;
  }

  return min;
}

template<typename T>
T getSum(vector<T> &vec)
{
  T sum = 0;
#pragma omp parallel for
  for (uint i=0; i<vec.size(); i++) {
    sum += vec[i];
  }

  return sum;
}

template<typename T>
void addVectors(vector<T> &vec1, vector<T> &vec2)
{
  if (vec1.size() != vec2.size()) FatalError("Vectors not of same size.");

  for (unsigned int i=0; i<vec1.size(); i++) vec1[i] += vec2[i];
}

template<typename T>
vector<T> & operator+=(vector<T>& lhs, vector<T>& rhs)
{
  if (lhs.size() != rhs.size()) FatalError("Vectors not of same size.");

  for (unsigned int i=0; i<lhs.size(); i++) lhs[i] += rhs[i];

  return lhs;
}

template<typename T>
vector<T> operator+(const vector<T>& lhs, vector<T>& rhs)
{
  if (lhs.size() != rhs.size()) FatalError("Vectors not of same size.");

  vector<T> out(lhs.size());
  for (unsigned int i=0; i<lhs.size(); i++) out[i] = lhs[i] + rhs[i];

  return out;
}

template<typename T>
vector<T> operator*(const vector<T>& lhs, double rhs)
{
  vector<T> out(lhs.size());
  for (unsigned int i=0; i<lhs.size(); i++) out[i] = lhs[i]*rhs;

  return out;
}

template<typename T>
vector<T> operator*(matrix<T>& mat, const vector<T> &vec)
{
  if (mat.getDim1() != vec.size()) FatalErrorST("Improper sizes for matrix*vector.");
  vector<T> out(mat.getDim0());
  for (int i=0; i<mat.getDim0(); i++) {
    for (int j=0; j<mat.getDim1(); j++) {
      out[i] += mat(i,j)*vec[j];
    }
  }
  return out;
}

Vec3 operator*(matrix<double>& mat, Vec3 &vec);

//----------Performance boost mod----------------------
/*template<typename T>
void std::vector<T>::operator*()*/
//----------Performance boost mod----------------------

template<typename T>
vector<T> operator/(const vector<T>& lhs, double rhs)
{
  vector<T> out(lhs.size());
  for (unsigned int i=0; i<lhs.size(); i++) out[i] = lhs[i]/rhs;

  return out;
}

class simTimer {
private:
  std::chrono::high_resolution_clock::time_point initTime;
  std::chrono::high_resolution_clock::time_point finalTime;

public:
  void startTimer();
  void stopTimer();
  void showTime(int precision=3);
  double getElapsedTime(void);
};

#ifdef _OMP
void omp_blocked_dgemm(CBLAS_ORDER mode, CBLAS_TRANSPOSE transA,
    CBLAS_TRANSPOSE transB, int M, int N, int K, double alpha, double* A, int lda,
    double* B, int ldb, double beta, double* C, int ldc);
#endif
