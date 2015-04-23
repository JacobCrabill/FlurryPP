/*!
 * \file global.cpp
 * \brief Definition of global constants, objects, and variables
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

#include "../include/global.hpp"

#include <cstdlib>
#include <string>

/* --- Misc. Common Constants --- */
double pi = 4.0*atan(1);

map<string,int> bcNum;

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void setGlobalVariables(void) {
  bcNum["none"] = NONE;
  bcNum["periodic"] = PERIODIC;
  bcNum["char"] = CHAR;
  bcNum["sup_in"] = SUP_IN;
  bcNum["sup_out"] = SUP_OUT;
  bcNum["sub_in"] = SUB_IN;
  bcNum["sub_out"] = SUB_OUT;
  bcNum["slip_wall"] = SLIP_WALL;
  bcNum["isothermal_noslip"] = ISOTHERMAL_NOSLIP;
  bcNum["adiabatic_noslip"] = ADIABATIC_NOSLIP;
}
