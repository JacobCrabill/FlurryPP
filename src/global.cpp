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

#ifndef _NO_MPI
#include "mpi.h"
#endif

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


void simTimer::startTimer(void)
{
  initTime = std::chrono::high_resolution_clock::now();
}

void simTimer::stopTimer(void)
{
  finalTime = std::chrono::high_resolution_clock::now();
}

void simTimer::showTime(void)
{
#ifndef _NO_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank == 0) {
#endif

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( finalTime - initTime ).count();
  double execTime = (double)duration/1000.;
  if (execTime > 60) {
    int minutes = floor(execTime/60);
    double seconds = execTime-(minutes*60);
    cout.setf(ios::fixed, ios::floatfield);
    cout << "Execution time = " << minutes << "min " << setprecision(3) << seconds << "s" << endl;
  }
  else {
    cout << setprecision(3) << "Execution time = " << execTime << "s" << endl;
  }

#ifndef _NO_MPI
  }
#endif
}

vector<int> getOrder(vector<double> &data)
{
  vector<pair<double,size_t> > vp;
  vp.reserve(data.size());
  for (size_t i = 0 ; i != data.size() ; i++) {
    vp.push_back(make_pair(data[i], i));
  }

  // Sorting will put lower values [vp.first] ahead of larger ones,
  // resolving ties using the original index [vp.second]
  sort(vp.begin(), vp.end());
  vector<int> ind(data.size());
  for (size_t i = 0 ; i != vp.size() ; i++) {
    ind[i] = vp[i].second;
  }

  return ind;
}

double getCFLLimit(int order)
{
  switch(order) {
    case 0:
      return 1.393;

    case 1:
      return 0.464;

    case 2:
      return 0.235;

    case 3:
      return 0.139;

    case 4:
      return 0.100;

    case 5:
      return 0.068;

    default:
      FatalError("CFL limit no available for this order!");
  }
}
