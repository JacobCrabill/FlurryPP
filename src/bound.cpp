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
  normL.setup(nFptsL,nDims);
  detJacL.resize(nFptsL);

  // Get access to data at left element
  int fpt=0;
  for (int i=fptStartL; i<fptEndL; i++) {
    UL[fpt] = (eL->U_fpts[i]);
    detJacL[fpt] = (eL->detJac_fpts[i]); // change to double**[] = &(eL->det[])

    for (int dim=0; dim<nDims; dim++)
      normL[fpt][dim] = (eL->norm_fpts[i][dim]); // change to dbl ptr

    FL[fpt].setup(nDims,nFields);
    for (int dim=0; dim<nDims; dim++)
      for (int k=0; k<nFields; k++)
        FL[fpt][dim][k] = &(eL->F_fpts[dim][i][k]);

    fpt++;
  }

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
  tempUL.resize(nFields);
}

void bound::calcInviscidFlux()
{
  // Set the boundary condition [store in UC]
  applyBCs();

  for (int i=0; i<nFptsL; i++) {
    for (int j=0; j<nFields; j++)
      tempUL[j] = UL[i][j]/(detJacL[i]);

    // Calcualte discontinuous inviscid flux at flux points
    inviscidFlux(tempUL,tempFL, params);

    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      upwindFlux(tempUL, UC[i], normL[i], Fn[i], params);
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
