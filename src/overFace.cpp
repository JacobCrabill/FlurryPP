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

      if (params->oversetMethod == 1)
        Fn(i,k) = OComm->U_in(fptOffset+i,nFields+k);
    }
  }
}

void overFace::getRightGradient(void)
{
  // Note: fptOffset must be set by Solver during overset setup
  for (int i=0; i<nFptsL; i++) {
    for (int dim=0; dim<nDims; dim++) {
      for (int k=0; k<nFields; k++) {
        gradUR[i](dim,k) = OComm->gradU_in(fptOffset+i,dim*nFields+k);
      }
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
  vector<double> force = {0,0,0,0,0,0};
  return force;
}

vector<double> overFace::computeMassFlux()
{
  // Not an inlet/outlet boundary - return 0
  vector<double> force(nFields);
  return force;
}

vector<point> overFace::getPosFpts()
{
  for (uint i=0; i<nFptsL; i++) {
    posFpts[i] = eL->getPosFpt(fptStartL+i);;
  }
  return posFpts;
}

vector<point> overFace::getNormFpts()
{
  vector<Vec3> normFpts(nFptsL);
  for (uint i=0; i<nFptsL; i++) {
    normFpts[i] = point(&eL->norm_fpts(fptStartL+i,0),nDims);
  }
  return normFpts;
}

void overFace::rusanovFlux(void)
{
  /* Different flux function to utilize interpolated corrected flux, but still
   * apply upwinding */
  if (params->oversetMethod == 1) {
    for (int fpt=0; fpt<nFptsL; fpt++) {
      // Get primitive variables
      double rhoL = UL(fpt,0);     double rhoR = UR(fpt,0);
      double uL = UL(fpt,1)/rhoL;  double uR = UR(fpt,1)/rhoR;
      double vL = UL(fpt,2)/rhoL;  double vR = UR(fpt,2)/rhoR;

      double wL, pL, vnL=0.;
      double wR, pR, vnR=0.;
      double vgn=0.;

      // Calculate pressure
      if (params->nDims==2) {
        pL = (params->gamma-1.0)*(UL(fpt,3)-0.5*rhoL*(uL*uL+vL*vL));
        pR = (params->gamma-1.0)*(UR(fpt,3)-0.5*rhoR*(uR*uR+vR*vR));
      }
      else {
        wL = UL(fpt,3)/rhoL;   wR = UR(fpt,3)/rhoR;
        pL = (params->gamma-1.0)*(UL(fpt,4)-0.5*rhoL*(uL*uL+vL*vL+wL*wL));
        pR = (params->gamma-1.0)*(UR(fpt,4)-0.5*rhoR*(uR*uR+vR*vR+wR*wR));
      }

      // Get normal fluxes, normal velocities
      for (int dim=0; dim<params->nDims; dim++) {
        vnL += normL(fpt,dim)*UL(fpt,dim+1)/rhoL;
        vnR += normL(fpt,dim)*UR(fpt,dim+1)/rhoR;
        if (params->motion)
          vgn += normL(fpt,dim)*Vg(fpt,dim);
      }

      // Get maximum eigenvalue for diffusion coefficient
      double csqL = max(params->gamma*pL/rhoL,0.0);
      double csqR = max(params->gamma*pR/rhoR,0.0);
      double eigL = std::fabs(vnL) + sqrt(csqL);
      double eigR = std::fabs(vnR) + sqrt(csqR);
      double eig  = max(eigL,eigR);

      /* Calculate Rusanov flux using corrected normal flux from donor grid
       * (copied to Fn) */

      // Outflow - use internal state
      if (vnL - vgn > 0) {
        inviscidFlux(UL[fpt],tempFL,params);
        double tempFnL[5] = {0,0,0,0,0};
        for (int k=0; k<nFields; k++)
          for (int dim=0; dim<nDims; dim++)
            tempFnL[k] += tempFL(dim,k)*normL(fpt,dim);

        for (int k=0; k<params->nFields; k++) {
          Fn(fpt,k) = tempFnL[k] - 0.5*eig*(UR(fpt,k)-UL(fpt,k));
        }
      }

      // Inflow - use external state [previously copied into Fn]
      else {
        for (int k=0; k<params->nFields; k++) {
          Fn(fpt,k) = Fn(fpt,k) - 0.5*eig*(UR(fpt,k)-UL(fpt,k));
        }
      }

      // Store wave speed for calculation of allowable dt
      if (params->motion) {
        eigL = std::fabs(vnL-vgn) + sqrt(csqL);
        eigR = std::fabs(vnR-vgn) + sqrt(csqR);
      }
      *waveSp[fpt] = max(eigL,eigR);
    }
  }

  // For other overset methods, just call normal Rusanov flux
  else {
    face::rusanovFlux();
  }
}
