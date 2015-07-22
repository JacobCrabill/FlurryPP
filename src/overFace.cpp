/*!
 * \file overFace.cpp
 * \brief Class to handle interface flux calculations on overset boundaries
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

#include "overFace.hpp"

#include "flux.hpp"
#include "ele.hpp"

void overFace::setupRightState(void)
{
  for (uint i=0; i<nFptsL; i++) {
    point pt = eL->getPosFpt(fptStartL+i);
    posFpts.push_back(pt);
  }

  // Get access to normal flux storage at right element [use look-up table to get right fpt]
  for (int i=0; i<nFptsL; i++) {
    FR[i] = (Solver->F_overPts[fptR[i]]);
  }
}

void overFace::getRightState(void)
{
  // Get data from right element [order reversed to match left ele]
  for (int fpt=0; fpt<nFptsL; fpt++) {
    for (int j=0; j<nFields; j++) {
      UR(fpt,j) = (eR->U_fpts(fptR[fpt],j));
    }

    // For dynamic grids, need to update geometry-related data
    if ((params->iter == params->initIter+1) || (params->motion != 0)) {
      for (int dim=0; dim<nDims; dim++) {
        normR(fpt,dim) = (eR->norm_fpts(fptR[fpt],dim));
      }
      dAR[fpt] = (eR->dA_fpts[fptR[fpt]]);
      detJacR[fpt] = (eR->detJac_fpts[fptR[fpt]]);
    }

    if (params->viscous) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUR[fpt](dim,j) = (eR->dU_fpts[dim](fptR[fpt],j));
    }
  }
}

void overFace::setRightStateFlux(void)
{
  for (int i=0; i<nFptsR; i++)
    for (int j=0; j<nFields; j++)
      FnR[i][j] = -Fn(i,j)*dAR[i]; // opposite normal direction
}

void overFace::setRightStateSolution(void)
{
  for (int i=0; i<nFptsR; i++)
    for (int j=0; j<nFields; j++)
      UcR[i][j] = UC(i,j);
}

vector<double> overFace::computeWallForce()
{
  // Not a wall boundary - return 0
  vector<double> force = {0,0,0};
  return force;
}

vector<point> overFace::getPosFpts()
{
  return posFpts;
}
