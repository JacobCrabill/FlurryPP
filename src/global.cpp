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

#include "global.hpp"

#include <cstdlib>
#include <string>

#ifndef _NO_MPI
#include "mpi.h"
#endif

#include "cblas.h"

/* --- Misc. Common Constants --- */
double pi = 4.0*atan(1);

//! Maps a boundary-condition string to its integer enum
// NOTE: 'symmetry' is just a psuedonym for 'slip_wall' which will not be
// considered a "wall" boundary condition for overset grids, force calc, etc.
map<string,int> bcStr2Num = {
  {"none", NONE},
  {"fluid", NONE},
  {"periodic", PERIODIC},
  {"char", CHAR},
  {"sup_in", SUP_IN},
  {"sup_out", SUP_OUT},
  {"sub_in", SUB_IN},
  {"sub_out", SUB_OUT},
  {"slip_wall", SLIP_WALL},
  {"isothermal_noslip", ISOTHERMAL_NOSLIP},
  {"adiabatic_noslip", ADIABATIC_NOSLIP},
  {"overset", OVERSET},
  {"symmetry", SYMMETRY}
};

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
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

std::ostream& operator<<(std::ostream &os, const point &pt)
{
  os << "(x,y,z) = " << pt.x << ", " << pt.y << ", " << pt.z;
  return os;
}

std::ostream& operator<<(std::ostream &os, const matrix<double> &mat)
{
  mat.print();
  return os;
}

std::ostream& operator<<(std::ostream &os, const vector<int> &vec)
{
  for (auto &val:vec) cout << val << ", ";
  return os;
}


std::ostream& operator<<(std::ostream &os, const vector<double> &vec)
{
  for (auto &val:vec) cout << val << ", ";
  return os;
}

std::ostream& operator<<(std::ostream &os, const matrix<int> &mat)
{
  mat.print();
  return os;
}

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

double simTimer::getElapsedTime(void)
{
  finalTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( finalTime - initTime ).count();
  return (double) duration / 1000.;
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

#ifdef _OMP
void omp_blocked_dgemm(CBLAS_ORDER mode, CBLAS_TRANSPOSE transA,
    CBLAS_TRANSPOSE transB, int M, int N, int K, double alpha, double* A, int lda,
    double* B, int ldb, double beta, double* C, int ldc)
{
  if (mode == CblasRowMajor)
  {
#pragma omp parallel
    {
      int nThreads = omp_get_num_threads();
      int thread_idx = omp_get_thread_num();

      int block_size = N / nThreads;
      int start_idx = block_size * thread_idx;

      if (thread_idx == nThreads-1)
        block_size += N % (block_size);

      cblas_dgemm(mode, transA, transB, M, block_size, K, alpha, A, lda,
                  B + start_idx, ldb, beta, C + start_idx, ldc);
    }
  }
  else
  {
#pragma omp parallel
    {
      int nThreads = omp_get_num_threads();
      int thread_idx = omp_get_thread_num();

      int block_size = N / nThreads;
      int start_idx = block_size * thread_idx;

      if (thread_idx == nThreads-1)
        block_size += N % (block_size);

      cblas_dgemm(mode, transA, transB, M, block_size, K, alpha, A, lda,
          B + ldb * start_idx, ldb, beta, C + ldc * start_idx, ldc);
    }
  }
}
#endif
