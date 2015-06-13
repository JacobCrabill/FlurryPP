/*!
 * \file face.cpp
 * \brief Class to handle interface flux calculations between elements
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

#include "../include/intFace.hpp"

#include "../include/flux.hpp"
#include "../include/ele.hpp"

void intFace::setupRightState(void)
{
  // This is kinda messy, but avoids separate initialize function
  locF_R = rightParam;

  nFptsR = eR->order+1;

  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartR = (locF_R*(nFptsR)) + nFptsR;
  fptEndR = (locF_R*(nFptsR));

  //FR.resize(nFptsR);
  FnR.resize(nFptsR);
  normR.setup(nFptsR,nDims);
  dAR.resize(nFptsR);
  detJacR.resize(nFptsL);

  if (params->viscous) {
    UcR.resize(nFptsR);
  }

  // Get access to normal flux storage at right element [order reversed to match left ele]
  int fpt = 0;
  for (int i=fptStartR-1; i>=fptEndR; i--) {
    FnR[fpt] = (eR->Fn_fpts[i]);

    if (params->viscous)
      UcR[fpt] = (eR->Uc_fpts[i]);

    fpt++;
  }
}

void intFace::getRightState(void)
{
  // Get data from right element [order reversed to match left ele]
  int fpt = 0;
  for (int i=fptStartR-1; i>=fptEndR; i--) {
    for (int j=0; j<nFields; j++) {
      UR(fpt,j) = (eR->U_fpts(i,j));
    }

    // For dynamic grids, need to update geometry-related data
    if ((params->iter == params->initIter+1) || (params->motion != 0)) {
      for (int dim=0; dim<nDims; dim++) {
        normR(fpt,dim) = (eR->norm_fpts(i,dim));
      }
      dAR[fpt] = (eR->dA_fpts[i]);
      detJacR[fpt] = (eR->detJac_fpts[i]);
    }

    if (params->viscous) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUR[fpt](dim,j) = (eR->dU_fpts[dim](i,j));
    }

    fpt++;
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
      UcR[i][j] = UC(i,j);
}
