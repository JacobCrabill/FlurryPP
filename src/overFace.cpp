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
  fptR.resize(nFptsL);

  for (uint i=0; i<nFptsL; i++) {
    point pt = eL->getPosFpt(fptStartL+i);
    posFpts.push_back(pt);
  }
}

void overFace::getRightState(void)
{
  // Note: fptOffset must be set by Solver during overset setup
  for (int i=0; i<nFptsL; i++) {
    for (int k=0; k<nFields; k++) {
      UR(i,k) = Solver->U_opts(fptOffset+i,k);
    }
  }
}

void overFace::setRightStateFlux(void)
{

}

void overFace::setRightStateSolution(void)
{

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
