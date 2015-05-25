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
#include "../include/ele.hpp"

void face::initialize(ele *eL, ele *eR, int locF_L, int rightParam, int gID, input *params)
{
  ID = gID;

  this->locF_L = locF_L;
  this->eL = eL;
  this->eR = eR;
  this->params = params;

  // Note: this is locF_R for internal faces, bcType for boundary faces
  this->rightParam = rightParam;

  nDims = params->nDims;
  nFields = params->nFields;

  // Setup temporary vectors for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);

  // Needed for MPI faces to separate communication from flux calculation
  isMPI = 0;
}

void face::setupFace(void)
{
  nFptsL = eL->order+1;

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.setup(nFptsL,nFields);
  UR.setup(nFptsL,nFields);
  FnL.resize(nFptsL);
  Fn.setup(nFptsL,nFields);
  normL.setup(nFptsL,nDims);
  dAL.resize(nFptsL);
  detJacL.resize(nFptsL);
  waveSp.resize(nFptsL);

  if (params->viscous) {
    // just a placeholder.  Need to properly size/reorder dimensions later.
    gradUL.resize(nDims);
    gradUR.resize(nDims);
    for (int i=0; i<nDims; i++) {
      gradUL[i].setup(nFptsL,nFields);
      gradUR[i].setup(nFptsR,nFields);
    }
  }

  // Get access to data at left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    FnL[fpt] = (eL->Fn_fpts[i]);
    waveSp[fpt] = &(eL->waveSp_fpts[i]);
    fpt++;
  }

  this->setupRightState();
}

void face::getLeftState()
{
  // Get data from left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    for (int j=0; j<nFields; j++) {
      UL(fpt,j) = (eL->U_fpts(i,j));
    }

    // For dynamic grids, need to update geometry-related data
    if ((params->iter == params->initIter+1) || (params->motion != 0)) {
      for (int dim=0; dim<nDims; dim++) {
        normL(fpt,dim) = (eL->norm_fpts(i,dim));
      }
      dAL[fpt] = (eL->dA_fpts[i]);
      detJacL[fpt] = (eL->detJac_fpts[i]);
    }

    fpt++;
  }
}

void face::calcInviscidFlux(void)
{
  if (!isMPI)
    getLeftState();  // Idea: instead of using ptrs to ele data, do copy?
  this->getRightState(); // <-- makes this more general for all face types, and allows face memory to be contiguous

  for (int i=0; i<nFptsL; i++) {
    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      laxFriedrichsFlux(UL[i], UR[i], normL[i], Fn[i], params);
    }
    else if (params->equation == NAVIER_STOKES) {
      if (params->riemann_type==0) {
        rusanovFlux(UL[i], UR[i], tempFL, tempFR, normL[i], Fn[i], waveSp[i], params);
      }
      else if (params->riemann_type==1) {
        roeFlux(UL[i], UR[i], normL[i], Fn[i], params);
      }
    }
  }

  // Transform normal flux using edge Jacobian and put into ele's memory
  for (int i=0; i<nFptsL; i++) {
    for (int j=0; j<nFields; j++)
      FnL[i][j] =  Fn(i,j)*dAL[i];

    *waveSp[i] /= dAL[i];
  }

  this->setRightState();
}

void face::calcViscousFlux(void)
{
  for (int i=0; i<nFptsL; i++) {
    // Calculate discontinuous viscous flux at flux points
    viscousFlux(UL[i], gradUL[i], tempFL, params);
    viscousFlux(UR[i], gradUR[i], tempFR, params);

    // Calculte common viscous flux at flux points
    ldgFlux(UL[i], UR[i], gradUL[i], gradUR[i], Fn[i], params);
  }
}
