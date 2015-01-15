/*!
 * \file bound.cpp
 * \brief Class to handle enforcement of boundary conditions at boundary faces
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

#include "../include/bound.hpp"


void bound::setupBound(ele *eL, int locF_L, int bcType, int gID)
{
  int fptStartL, fptEndL;

  ID = gID;
  this->bcType = bcType;
  this->locF_L = locF_L;

  nDims = params->nDims;
  nFields = params->nFields;

  nFptsL = eL->order+1;

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.resize(nFptsL);
  FL.resize(nFptsL);

  // Get access to data at left element
  for (int i=fptStartL; i<fptEndL; i++) {
    UL[i] = &(eL->U_fpts[i]);
    FL[i] = &(eL->F_fpts[i]);
    normL[i] = (eL->norm_fpts[i]);
  }

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
}

void bound::calcInviscidFlux()
{

}

void bound::calcViscousFlux()
{

}
