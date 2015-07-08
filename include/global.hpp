/*!
 * \file global.hpp
 * \brief Header file for global constants, objects, and variables
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

#include <omp.h>

#include "error.hpp"
#include "matrix.hpp"

// Forward declarations of basic Flurry classes
class geo;
class ele;
class solver;

using namespace std;

typedef unsigned int uint;

/* --- Misc. Common Constants / Globally-Useful Variables --- */
extern double pi;

extern map<string,int> bcStr2Num; //! Maps a boundary-condition string to its integer enum

/*! enumeration for element type */
enum ETYPE {
  TRI     = 0,
  QUAD    = 1,
  TET     = 2,
  PRISM   = 3,
  HEX     = 4,
  PYRAMID = 5
};

/*! Enumeration for mesh (either create cartesian mesh or read from file) */
enum MESH_TYPE {
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
  CHAR = 1,
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
  SYMMETRY = 12
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

  point(double* pt) {
    x = pt[0];
    y = pt[1];
    z = pt[2];
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

  double operator*(point b) {
    return x*b.x + y*b.y + z*b.z;
  }

  void abs(void) {
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);
  }

};

//! For clearer notation when a vector is implied, rather than a point
typedef struct point Vec3;

int factorial(int n);

void setGlobalVariables(void);

bool checkNaN(vector<double> &vec);

bool checkNaN(double* vec, int size);

/*! Get polynomial-order-based CFL limit.  Borrowed from Josh's zefr code. */
double getCFLLimit(int order);

//double randRange(double xMin, double xMax);

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

// Good for numeric types
template<typename T>
T getMax(vector<T> &vec)
{
  T max = 0;
  for (auto& i:vec) {
    if (i>max) max = i;
  }

  return max;
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
  void showTime();

};
