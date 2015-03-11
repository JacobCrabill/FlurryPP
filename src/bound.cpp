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
  // Set the boundary condition [store in UC]
  applyBCs();

  for (int i=0; i<nFptsL; i++) {


    // Calcualte discontinuous inviscid flux at flux points
    inviscidFlux(*UL[i],tempFL, params);

    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      upwindFlux(*UL[i], UC[i], normL[i], *Fn[i], params);
    }
//    else if (params->equation == NAVIER_STOKES) {
//      if (params->riemann_type==0) {
//        rusanovFlux(*UL[i], *UC[i], *FL[i], *FR[i], normL[i], *Fn[i], params);
//      }
//      else if (params->riemann_type==1) {
//        roeFlux(*UL[i], *UR[i], normL[i], *Fn[i], params);
//      }
//    }
  }
}

void bound::calcViscousFlux()
{

}

void bound::applyBCs(void)
{
  for (int fpt=0; fpt<nFptsL; fpt++) {
    if (params->equation == NAVIER_STOKES) {
      if (bcType == 0) {

      }
    }
    else if (params->equation == ADVECTION_DIFFUSION) {

    }
  }
}
