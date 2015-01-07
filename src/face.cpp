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

#include "../include/face.hpp"
#include "../include/flux.hpp"

face::face()
{

}

void face::setupFace(ele *eL, ele *eR, int locF_L, int locF_R, int gID)
{
  int fptStartL, fptEndL, fptStartR, fptEndR, i;

  ID = gID;

  nDims = params->nDims;
  nFields = params->nFields;

  nFptsL = eL->order+1;
  nFptsR = eR->order+1;

  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;
  fptStartR = (locF_R*(nFptsR)) + nFptsR;
  fptEndR = (locF_R*(nFptsR));

  UL.resize(nFptsL);
  UR.resize(nFptsR);
  FL.resize(nFptsL);
  FR.resize(nFptsR);

  // Get access to data at left element
  for (i=fptStartL; i<fptEndL; i++) {
    UL[i] = &(eL->U_fpts[i]);
    FL[i] = &(eL->F_fpts[i]);
    normL[i] = &(eL->norm_fpts[i]);
  }

  // Get access to data at right element [order reversed to match left ele]
  for (i=fptStartR-1; i>=fptEndR; i--) {
    UR[i] = &(eR->U_fpts[i]);
    FR[i] = &(eR->F_fpts[i]);
    normR[i] = &(eR->norm_fpts[i]);
  }

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
}

void face::calcInviscidFlux(void)
{
  int i, j, k;

  for (i=0; i<nFptsL; i++) {
    // Calcualte discontinuous inviscid flux at flux points
    inviscidFlux(*UL[i],tempFL, params);
    inviscidFlux(*UR[i],tempFR, params);

    // Calculate common inviscid flux at flux points
    if (params->riemann_type==0) {
      rusanovFlux(*UL[i], *UR[i], *FL[i], *FR[i], *normL[i], *Fn[i], params);
    }else if (params->riemann_type==1) {
      roeFlux(*UL[i], *UR[i], *normL[i], *Fn[i], params);
    }
  }
}

void face::calcViscousFlux(void)
{
  int i;

  for (i=0; i<nFptsL; i++) {
    // Calculate discontinuous viscous flux at flux points
    viscousFlux(*UL[i], *gradUL[i], tempFL, params);
    viscousFlux(*UR[i], *gradUR[i], tempFR, params);

    // Calculte common viscous flux at flux points
    ldgFlux(*UL[i], *UR[i], *gradUL[i], *gradUR[i], *Fn[i], params);
  }
}
