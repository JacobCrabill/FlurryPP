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
    //FR[i] = (Solver->F_opts[fptR[i]]);
  }
}

void overFace::getRightState(void)
{

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
