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
#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>    // std::size_t
#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>

#include "matrix.hpp"

using namespace std;

//! Prints the error message, the stack trace, and exits
#define FatalError(s) {                                             \
  printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        \
  exit(1); }


/* --- Misc. Common Constants --- */
extern double pi;

/** enumeration for element type */
enum ETYPE {
  TRI     = 0,
  QUAD    = 1,
  TET     = 2,
  PRISM   = 3,
  HEX     = 4,
  PYRAMID = 5
};

enum MESH_TYPE {
  READ_MESH   = 0,
  CREATE_MESH = 1
};

struct point
{
  double x, y, z;
};

/*! Find indices of all values in vec equal to val */
template<typename T>
vector<int> findEq(const vector<T> &vec, T val)
{
  vector<int> out;

  for (unsigned int i=0; i<vec.size(); i++) {
    if (vec[i]==val) out.push_back(i);
  }

  return out;
}

/*! Find index of first occurance of val in vec */
template<typename T>
int findFirst(const vector<int> &vec, T val)
{
  if (vec.size()==0) return -1;

  for (int i=0; i<(int)vec.size(); i++) {
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

