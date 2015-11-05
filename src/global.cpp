/*!
 * \file global.cpp
 * \brief Definition of global constants, objects, and variables
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "../include/global.hpp"

#include <cstdlib>
#include <string>

#ifndef _NO_MPI
#include "mpi.h"
#endif

/* --- Misc. Common Constants --- */
double pi = 4.0*atan(1);

map<string,int> bcStr2Num;

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void setGlobalVariables(void) {
  bcStr2Num["none"] = NONE;
  bcStr2Num["fluid"] = NONE;
  bcStr2Num["periodic"] = PERIODIC;
  bcStr2Num["char"] = CHAR;
  bcStr2Num["sup_in"] = SUP_IN;
  bcStr2Num["sup_out"] = SUP_OUT;
  bcStr2Num["sub_in"] = SUB_IN;
  bcStr2Num["sub_out"] = SUB_OUT;
  bcStr2Num["slip_wall"] = SLIP_WALL;
  bcStr2Num["isothermal_noslip"] = ISOTHERMAL_NOSLIP;
  bcStr2Num["adiabatic_noslip"] = ADIABATIC_NOSLIP;
  bcStr2Num["overset"] = OVERSET;
  bcStr2Num["symmetry"] = SYMMETRY;
  // NOTE: 'symmetry' is just a psuedonym for 'slip_wall' which will not be
  // considered a "wall" boundary condition for overset grids, force calc, etc.
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
      FatalError("CFL limit not available for this order!");
  }
}

vector<double> getQuadratureWeights1D(int order)
{
  // Order here refers to the order of a polynomial fit through
  // the Gauss points, not the order of accuracy of integration
  // using the same number of points

  vector<double> outWts(order+1);

  if (order == 0) {
    outWts[0] =  2.0;
  }
  else if(order == 1) {
    outWts[0] = 1.0;
    outWts[1] = 1.0;
  }
  else if(order == 2) {
    outWts[0] = 0.5555555555555556;
    outWts[1] = 0.8888888888888888;
    outWts[2] = 0.5555555555555556;
  }
  else if(order == 3) {
    outWts[0] = 0.3478548451374538;
    outWts[1] = 0.6521451548625461;
    outWts[2] = 0.6521451548625461;
    outWts[3] = 0.3478548451374538;
  }
  else if(order == 4) {
    outWts[0] = 0.2369268850561891;
    outWts[1] = 0.4786286704993665;
    outWts[2] = 0.000000000000000;
    outWts[3] = 0.4786286704993665;
    outWts[4] = 0.2369268850561891;
  }
  else {
    cout << "Order = " << order << ": " << flush;
    FatalError("Gauss quadrature weights for this order not implemented.");
  }

  return outWts;
}

Vec3 operator*(matrix<double>& mat, Vec3 &vec)
{
  Vec3 out;
  int nDims = mat.getDim1();
  for (int i=0; i<nDims; i++) {
    for (int j=0; j<nDims; j++) {
      out[i] += mat(i,j)*vec[j];
    }
  }
  return out;
}

point operator/(point a, double b) { return a/=b; }
point operator*(point a, double b) { return a*=b; }


bool operator<(const point &a, const point &b) { return a.x < b.x; }

double getDist(point a, point b)
{
  Vec3 dx = a - b;
  return dx.norm();
}

matrix<double> createMatrix(vector<point> &pts)
{
  matrix<double> out(pts.size(),3);

  for (uint i=0; i<pts.size(); i++) {
    for (uint j=0; j<3; j++) {
      out(i,j) = pts[i][j];
    }
  }

  return out;
}

void simTimer::startTimer(void)
{
  initTime = std::chrono::high_resolution_clock::now();
}

void simTimer::stopTimer(void)
{
  finalTime = std::chrono::high_resolution_clock::now();
}

void simTimer::showTime(int precision)
{
  int rank = 0;
#ifndef _NO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  if (rank == 0) {
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( finalTime - initTime ).count();
    double execTime = (double)duration/1000.;
    cout.setf(ios::fixed, ios::floatfield);
    if (execTime > 60) {
      int minutes = floor(execTime/60);
      double seconds = execTime-(minutes*60);
      cout << "Execution time = " << minutes << "min " << setprecision(precision) << seconds << "s" << endl;
    }
    else {
      cout << setprecision(precision) << "Execution time = " << execTime << "s" << endl;
    }
  }
}
