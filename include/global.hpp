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
#include <vector>
#include <stdio.h>

//#include "matrix.hpp"

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

struct point
{
  double x, y, z;
};

template <typename T>
typedef vector<vector<double>> matrix;
