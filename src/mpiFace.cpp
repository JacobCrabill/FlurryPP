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

#include "../include/mpiFace.hpp"

#ifndef _NO_MPI
#include <mpi.h>
#endif

#include "../include/flux.hpp"
#include "../include/ele.hpp"

void mpiFace::setupRightState(void)
{
#ifdef _NO_MPI
  FatalError("Trying to setup an MPI face but code is not compiled for MPI.");
#endif

  // TODO: set procR, idR, etc.
  procR = rightParam;


  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartR = (locF_R*(nFptsR)) + nFptsR;
  fptEndR = (locF_R*(nFptsR));

  UR.setup(nFptsR,nFields);
  FR.resize(nFptsR);
  FnR.resize(nFptsR);
  normR.setup(nFptsR,nDims);
  dAR.resize(nFptsR);
  detJacR.resize(nFptsL);

  // Get access to normal flux storage at right element [order reversed to match left ele]
  int fpt = 0;
  for (int i=fptStartR-1; i>=fptEndR; i--) {
    FnR[fpt] = (eR->Fn_fpts[i]);
    FR[fpt].setup(nDims,nFields);
    fpt++;
  }
}

void mpiFace::getRightState(void)
{
  // Get data from right element [order reversed to match left ele]

}

void mpiFace::setRightState(void)
{
  /* Options:
   * 1) Have duplicate mpiFaces (one on either side of mpi boundary);
   *    each one sends left state to the other's right state here
   * 2) Don't duplicate mpiFaces; have it set FnR in the right ele [but how
   *    / where to put MPI calls?]
   */
  // -- Going with Option 1 --
  // Create outgoing, incoming buffers

}
