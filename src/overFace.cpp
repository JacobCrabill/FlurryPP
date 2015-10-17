/*!
 * \file overFace.cpp
 * \brief Class to handle interface flux calculations on overset boundaries
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */

#include "overFace.hpp"

#include "flux.hpp"
#include "ele.hpp"
#include "overComm.hpp"

void overFace::setupRightState(void)
{
  for (uint i=0; i<nFptsL; i++) {
    point pt = eL->getPosFpt(fptStartL+i);
    posFpts.push_back(pt);
  }
}

void overFace::getPointersRight()
{
  // Do nothing
}

void overFace::getRightState(void)
{
  // Note: fptOffset must be set by Solver during overset setup
  for (int i=0; i<nFptsL; i++) {
    for (int k=0; k<nFields; k++) {
      UR(i,k) = OComm->U_in(fptOffset+i,k);
    }
  }
}

void overFace::getRightGradient(void)
{
  // Note: fptOffset must be set by Solver during overset setup
//  for (int i=0; i<nFptsL; i++) {
//    for (int dim=0; dim<nDims; dim++) {
//      for (int k=0; k<nFields; k++) {
//        gradUR(i,dim,k) = OComm->gradU_in(fptOffset+i,dim,k);
//      }
//    }
//  }
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
  vector<double> force = {0,0,0,0,0,0};
  return force;
}

vector<point> overFace::getPosFpts()
{
  for (uint i=0; i<nFptsL; i++) {
    posFpts[i] = eL->getPosFpt(fptStartL+i);;
  }
  return posFpts;
}
