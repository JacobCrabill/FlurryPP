/*!
 * \file face.cpp
 * \brief Class to handle interface flux calculations between elements
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

#include "intFace.hpp"

#include "flux.hpp"
#include "ele.hpp"

void intFace::setupRightState(void)
{
  // This is kinda messy, but avoids separate initialize function
  faceID_R = myInfo.IDR;
  relRot = myInfo.relRot;

  if (nDims == 2)
    nFptsR = eR->order+1;
  else
    nFptsR = (eR->order+1)*(eR->order+1);

  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  FnR.resize(nFptsR);
  normR.setup(nFptsR,nDims);
  dAR.resize(nFptsR);

  /* --- Setup the L/R flux-point matching --- */
  fptR.resize(nFptsL);
  if (nDims == 2) {
    // For 1D faces [line segments] only - find first/last ID of fpts;
    // right element's are simply reversed
    fptStartR = (faceID_R*(nFptsR)) + nFptsR;
    fptEndR = (faceID_R*(nFptsR));

    int fpt = 0;
    for (int i=fptStartR-1; i>=fptEndR; i--) {
      fptR[fpt] = i;
      fpt++;
    }
  }
  else if (nDims == 3) {
    // Only for quad tensor-product faces
    int order = sqrt(nFptsL)-1;

    fptStartR = nFptsR*faceID_R;
    fptEndR = nFptsR*(faceID_R+1);

    for (int i=0; i<nFptsL; i++) {
      int ifpt = i%(order+1);
      int jfpt = floor(i/(order+1));
      switch (relRot) {
        case 0:
          fptR[i] = fptStartR + ifpt*(order+1) + jfpt;
          break;
        case 1:
          fptR[i] = fptStartR + order - ifpt + jfpt*(order+1);
          break;
        case 2:
          fptR[i] = fptEndR - 1 - (ifpt*(order+1) + jfpt);
          break;
        case 3:
          fptR[i] = fptEndR - (order+1)*(jfpt+1) + ifpt;
          break;
      }
    }
  }

  if (params->viscous) {
    dUcR.setup(nFptsR,nFields);
  }
}

void intFace::getPointersRight(void)
{
  // Get access to normal flux storage at right element [use look-up table to get right fpt]
  for (int i=0; i<nFptsL; i++) {
    FnR[i] = &(eR->Fn_fpts(fptR[i],0));

    if (params->viscous)
      for (int k = 0; k < nFields; k++)
        dUcR(i,k) = &(eR->dUc_fpts(fptR[i],k));
  }
}

void intFace::getRightState(void)
{
  // Get data from right element [order reversed to match left ele]
  for (int fpt=0; fpt<nFptsL; fpt++) {
    for (int j=0; j<nFields; j++) {
      UR(fpt,j) = (eR->U_fpts(fptR[fpt],j));
    }

    /* For dynamic grids (besides rigid translation), need to update
     * geometry-related data on every iteration, not just during setup */
    if (isNew_R || (params->motion != 0 && params->motion != 4)) {
      for (int dim=0; dim<nDims; dim++) {
        normR(fpt,dim) = (eR->norm_fpts(fptR[fpt],dim));
      }
      dAR[fpt] = (eR->dA_fpts(fptR[fpt]));
    }
  }

  isNew_R = false;
}

void intFace::getRightGradient(void)
{
  // Get data from right element [order reversed to match left ele]
  if (params->viscous) {
    for (int fpt=0; fpt<nFptsL; fpt++) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUR[fpt](dim,j) = (eR->dU_fpts(dim,fptR[fpt],j));
    }
  }
}

void intFace::setRightStateFlux(void)
{
  for (int i=0; i<nFptsR; i++)
    for (int j=0; j<nFields; j++)
      FnR[i][j] = -Fn(i,j)*dAR[i]; // opposite normal direction


}

void intFace::setRightStateSolution(void)
{
  for (int i=0; i<nFptsR; i++)
    for (int j=0; j<nFields; j++)
      *dUcR(i,j) = UC(i,j) - UR(i,j);
}

vector<double> intFace::computeWallForce()
{
  // Not a wall boundary - return 0
  vector<double> force = {0,0,0,0,0,0};
  return force;
}

vector<double> intFace::computeMassFlux()
{
  // Not an inlet/outlet boundary - return 0
  vector<double> force(nFields);
  return force;
}


void intFace::get_U_index(int fpt, int& ind, int& stride)
{
  /* U : nFpts x nVars x 2 */
  int ic1 = eL->ID;
  int ic2 = eR->ID;

  if (ic1 > 0 && Solver->Geo->iblankCell[ic1] == NORMAL)
  {
    ind    = std::distance(&Solver->U_fpts(0,0,0), &eR->U_fpts(fptStartL+fpt,0));
    stride = std::distance(&eR->U_fpts(fptStartL+fpt,0), &eR->U_fpts(fptStartL+fpt,1));
  }
  else if (ic2 > 0 && Solver->Geo->iblankCell[ic2] == NORMAL)
  {
    ind    = std::distance(&Solver->U_fpts(0,0,0), &eL->U_fpts(fptStartL+fpt,0));
    stride = std::distance(&eL->U_fpts(fptStartL+fpt,0), &eL->U_fpts(fptStartL+fpt,1));
  }
  else
  {
    //printf("face %d: ibf %d | ic1,2: %d,%d, ibc1,2: %d,%d\n",faceID,geo->iblank_face[faceID],ic1,ic2,geo->iblank_cell[ic1],geo->iblank_cell[ic2]);
    FatalError("get_U_index : Face not blanked but both elements are!");
  }
}
