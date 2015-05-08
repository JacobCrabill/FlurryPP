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

void face::setupFace(ele *eL, ele *eR, int locF_L, int locF_R, int gID, input *params)
{
  int fptStartL, fptEndL, fptStartR, fptEndR, i, fpt;

  ID = gID;

  this->locF_L = locF_L;
  this->locF_R = locF_R;
  this->eL = eL;
  this->eR = eR;
  this->params = params;

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
  FnL.resize(nFptsL);
  FnR.resize(nFptsR);
  Fn.setup(nFptsL,nFields);
  normL.resize(nFptsL);//,nDims);
  normR.resize(nFptsR);//,nDims);
  dAL.resize(nFptsL);
  dAR.resize(nFptsR);
  detJacL.resize(nFptsL);
  detJacR.resize(nFptsL);

  posFpts.resize(nFptsL); // Probably only needed for debugging.  Remove later.

  // Get access to data at left element
  fpt = 0;
  for (i=fptStartL; i<fptEndL; i++) {
    UL[fpt] = (eL->U_fpts[i]);
    FnL[fpt] = (eL->Fn_fpts[i]);
    dAL[fpt] = &(eL->dA_fpts[i]);
    detJacL[fpt] = (eL->detJac_fpts[i]);
    posFpts[fpt] = eL->pos_fpts[i];
    normL[fpt] = (eL->norm_fpts[i]);

    FL[fpt].setup(nDims,nFields);
    for (int dim=0; dim<nDims; dim++)
      for (int k=0; k<nFields; k++)
        FL[fpt][dim][k] = &(eL->F_fpts[dim][i][k]);

    fpt++;
  }

  // Get access to data at right element [order reversed to match left ele]
  fpt = 0;
  for (i=fptStartR-1; i>=fptEndR; i--) {
    UR[fpt] = (eR->U_fpts[i]);
    FnR[fpt] = (eR->Fn_fpts[i]);
    dAR[fpt] = &(eR->dA_fpts[i]);
    detJacR[fpt] = (eR->detJac_fpts[i]);
    normR[fpt] = (eR->norm_fpts[i]);

    FR[fpt].setup(nDims,nFields);
    for (int dim=0; dim<nDims; dim++)
      for (int k=0; k<nFields; k++)
        FR[fpt][dim][k] = &(eR->F_fpts[dim][i][k]);

    fpt++;
  }

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);
  tempUR.resize(nFields);
}

void face::calcInviscidFlux(void)
{
  for (int i=0; i<nFptsL; i++) {
    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      laxFriedrichsFlux(UL[i], UR[i], normL[i], Fn[i], params);
    }
    else if (params->equation == NAVIER_STOKES) {
      if (params->riemann_type==0) {
        // Calcualte discontinuous inviscid flux at flux points
        inviscidFlux(UL[i], tempFL, params);
        inviscidFlux(UR[i], tempFR, params);
        rusanovFlux(UL[i], UR[i], tempFL, tempFR, normL[i], Fn[i], params);
      }
      else if (params->riemann_type==1) {
        roeFlux(UL[i], UR[i], normL[i], Fn[i], params);
      }
    }

    // Calculate difference between discontinuous & common normal flux, and store in ele
    // (Each ele needs only the difference, not the actual common value, for the correction)
    // Need dAL/R to transform normal flux back to reference space
    for (int j=0; j<nFields; j++) {
      FnL[i][j] =  Fn[i][j]*(*dAL[i]);
      FnR[i][j] = -Fn[i][j]*(*dAR[i]); // opposite normal direction
    }
  }
}

void face::calcViscousFlux(void)
{
  int i;

  for (i=0; i<nFptsL; i++) {
    // Calculate discontinuous viscous flux at flux points
    viscousFlux(UL[i], *gradUL[i], tempFL, params);
    viscousFlux(UR[i], *gradUR[i], tempFR, params);

    // Calculte common viscous flux at flux points
    ldgFlux(UL[i], UR[i], *gradUL[i], *gradUR[i], Fn[i], params);
  }
}
