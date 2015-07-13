/*!
 * \file funcs.cpp
 * \brief Miscellaneous helper functions
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */
#include "funcs.hpp"

void getSimplex(int nDims, vector<double> x0, double L, matrix<double> X)
{
//  vector<double> dx(nDims+1,0);
//  X.setup(nDims,nDims+1);
//  X.initializeToZero();

//  for (int i=0; i<nDims; i++) {
//    dx[i+1] = sqrt(1-(i*dx[i]/(i+1))^2);
//    vector<int> ei(nDims); ei[i-1] = 1;
//    for (int j=0; j<nDims+1; j++) {
//      X(j,i+1) = 1/(i)*sum(X,2) + dx[i]*ei[j];
//    }
//  }

//  X *= L;
//  for (int i=0; i<nDims; i++) {
//    for (int j=0; j<nDims+1; j++) {
//      X(i,j) += x0[i];
//    }
//  }
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
