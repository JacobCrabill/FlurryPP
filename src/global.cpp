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
#include <cstring>
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
  bcNum["fluid"] = NONE;
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


bool checkNaN(vector<double> &vec)
{
  for (auto& i:vec)
    if (std::isnan(i)) return true;

  return false;
}

bool checkNaN(double* vec, int size)
{
  for (int i=0; i<size; i++)
    if (std::isnan(vec[i])) return true;

  return false;
}

void bitswap(int value, char *buf)
{
  union temp
  {
    int  value;
    char c[4];
  } in, out;
  in.value = value;
  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  memcpy(buf, out.c, 4);
}

void bitswap(uint value, char *buf)
{
  union temp
  {
    uint  value;
    char c[4];
  } in, out;
  in.value = value;
  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  memcpy(buf, out.c, 4);
}

void bitswap(float value, char *buf)
{
  union temp
  {
    float value;
    char  c[4];
  } in, out;
  in.value = value;
  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  memcpy(buf, out.c, 4);
}

void bitswap(double value, char *buf)
{
  union temp
  {
    double value;
    char   c[8];
  } in, out;
  in.value = value;
  out.c[0] = in.c[7];
  out.c[1] = in.c[6];
  out.c[2] = in.c[5];
  out.c[3] = in.c[4];
  out.c[4] = in.c[3];
  out.c[5] = in.c[2];
  out.c[6] = in.c[1];
  out.c[7] = in.c[0];
  memcpy(buf, out.c, 8);
}

void bitunswap(char *buf, int *value)
{
  union temp
  {
    int  value;
    char c[4];
  } in, out;
  memcpy(in.c, buf, 4);
  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  memcpy(value, &out.value, 4);
}

void bitunswap(char *buf, float *value)
{
  union temp
  {
    float value;
    char  c[4];
  } in, out;
  memcpy(in.c, buf, 4);
  out.c[0] = in.c[3];
  out.c[1] = in.c[2];
  out.c[2] = in.c[1];
  out.c[3] = in.c[0];
  memcpy(value, &out.value, 4);
}

void bitunswap(char *buf, double *value)
{
  union temp
  {
    double value;
    char   c[8];
  } in, out;
  memcpy(in.c, buf, 8);
  out.c[0] = in.c[7];
  out.c[1] = in.c[6];
  out.c[2] = in.c[5];
  out.c[3] = in.c[4];
  out.c[4] = in.c[3];
  out.c[5] = in.c[2];
  out.c[6] = in.c[1];
  out.c[7] = in.c[0];
  memcpy(value, &out.value, 8);
}
