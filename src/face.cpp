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
  int fptStartL, fptEndL, fptStartR, fptEndR, i, fpt;

  ID = gID;

  this->locF_L = locF_L;
  this->locF_R = locF_R;

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
  dFnL.resize(nFptsL);
  dFnR.resize(nFptsR);
  Fn.setup(nFptsL,nFields);
  normL.setup(nFptsL,nDims);
  normR.setup(nFptsR,nDims);
  dAL.resize(nFptsL);
  dAR.resize(nFptsR);
  detJacL.resize(nFptsL);
  detJacR.resize(nFptsL);

  posFpts.resize(nFptsL); // Probably only needed for debugging.  Remove later.

  // Get access to data at left element
  fpt = 0;
  for (i=fptStartL; i<fptEndL; i++) {
    /* Old Versions of FL, FR */
    //FR[fpt] = &(eR->F_fpts[i]);
    //(*FL[fpt]).setup(nDims,nFields);
    //for (int j=0; j<nDims; j++)
    //   FL[j][fpt] = &(eL->F_fpts[j][i]);
    UL[fpt] = (eL->U_fpts[i]);

    dFnL[fpt] = (eL->dFn_fpts[i]);
    dAL[fpt] = (eL->dA_fpts[i]);
    detJacL[fpt] = (eL->detJac_fpts[i]);
    posFpts[fpt] = eL->pos_fpts[i];

    for (int dim=0; dim<nDims; dim++)
      normL[fpt][dim] = (eL->norm_fpts[i][dim]);

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
    dFnR[fpt] = (eR->dFn_fpts[i]);
    dAR[fpt] = (eR->dA_fpts[i]);
    detJacR[fpt] = (eR->detJac_fpts[i]);

    for (int dim=0; dim<nDims; dim++)
      normR[fpt][dim] = (eR->norm_fpts[i][dim]);  // change norm to matrix<double*> for future

    FR[fpt].setup(nDims,nFields);
    for (int dim=0; dim<nDims; dim++)
      for (int k=0; k<nFields; k++)
        FR[fpt][dim][k] = &(eR->F_fpts[dim][i][k]);

    fpt++;
  }

  /* Other FL,FR storage methods
  // Trying out new method of storing FL, FR...
//  FL.resize(nFptsL);
//  FR.resize(nFptsR);
//  fpt = 0;
//  for (i=fptStartL; i<fptEndL; i++) {
//    FL[fpt].resize(nDims);
//    for (int j=0; j<nDims; j++) {
//      FL[fpt][j].resize(nFields);
//      for (int k=0; k<nFields; k++) {
//        FL[fpt][j][k] = &(eL->F_fpts[j][i][k]);
//      }
//    }
//    fpt++;
//  }
//
//  fpt = 0;
//  for (i=fptStartR-1; i>=fptEndR; i--) {
//    FR[fpt].resize(nDims);
//    for (int j=0; j<nDims; j++) {
//      FR[fpt][j].resize(nFields);
//      for (int k=0; k<nFields; k++) {
//        FR[fpt][j][k] = &(eR->F_fpts[j][i][k]);
//      }
//    }
//    fpt++;
//  }

  // other method...
//      FL[j][fpt] = &(eL->F_fpts[j][i]);
//  fpt = 0;
//  for (int j=0; j<nDims; j++) {
//    FR[j].resize(nDims);
//    for (i=fptStartR-1; i>=fptEndR; i--) {
//      //FR[j][fpt] = &(eR->F_fpts[j][i]);
//      FR[j][fpt].resize(nFields);
//      for (int k=0; k<nFields; k++) {
//        FR[j][fpt][k] = &(eR->F_fpts[j][i][k]);
//      }
//    }
//    fpt++;
//  }
*/

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);
  tempUR.resize(nFields);
}

void face::calcInviscidFlux(void)
{
  int i,j,k;
  double tempFnL, tempFnR;

  for (i=0; i<nFptsL; i++) {
    for (j=0; j<nFields; j++) {
      tempUL[j] = UL[i][j]/(detJacL[i]);
      tempUR[j] = UR[i][j]/(detJacR[i]);
    }

    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      laxFriedrichsFlux(tempUL, tempUR, normL[i], Fn[i], params);
    }
    else if (params->equation == NAVIER_STOKES) {
      if (params->riemann_type==0) {
        // Calcualte discontinuous inviscid flux at flux points
        inviscidFlux(tempUL, tempFL, params);
        inviscidFlux(tempUR, tempFR, params);
        //rusanovFlux(tempUL, tempUR, *FL[i], *FR[i], normL[i], Fn[i], params);
        rusanovFlux(tempUL, tempUR, tempFL, tempFR, normL[i], Fn[i], params);
      }
      else if (params->riemann_type==1) {
        roeFlux(tempUL, tempUR, normL[i], Fn[i], params);
      }
    }

    // Calculate difference between discontinuous & common normal flux, and store in ele
    // (Each ele needs only the difference, not the actual common value, for the correction)
    for (j=0; j<nFields; j++) {
      tempFnL =  Fn[i][j];
      tempFnR = -Fn[i][j]; // opposite normal direction
      for (k=0; k<nDims; k++) {
        tempFnL -= *(FL[i][k][j])*normL[i][k];
        tempFnR -= *(FR[i][k][j])*normR[i][k];
      }
      // Transform back to reference space & store in element
      dFnL[i][j] = tempFnL*dAL[i];
      dFnR[i][j] = tempFnR*dAR[i];
    }
  }
}

void face::calcViscousFlux(void)
{
  int i;

  for (i=0; i<nFptsL; i++) {
    for (int k=0; k<nFields; k++) {
      tempUL[k] = UL[i][k]/(detJacL[i]);
      tempUR[k] = UR[i][k]/(detJacR[i]);
    }

    // Calculate discontinuous viscous flux at flux points
    viscousFlux(tempUL, *gradUL[i], tempFL, params);
    viscousFlux(tempUR, *gradUR[i], tempFR, params);

    // Calculte common viscous flux at flux points
    ldgFlux(tempUL, tempUR, *gradUL[i], *gradUR[i], Fn[i], params);
  }
}
