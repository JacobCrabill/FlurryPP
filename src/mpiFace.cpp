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
#include "mpi.h"
#endif

#include "../include/flux.hpp"
#include "../include/ele.hpp"

void mpiFace::setupRightState(void)
{
#ifdef _NO_MPI
  FatalError("Trying to setup an MPI face but code is not compiled for MPI.");
#endif

#ifndef _NO_MPI
  IDR = rightParam;

  // Send/Get # of flux points to/from right element [order reversed to match left ele]
  MPI_Isend(&nFptsL,1,MPI_INT,procR,IDR,MPI_COMM_WORLD,&nFpts_out);
  MPI_Irecv(&nFptsR,1,MPI_INT,procR,ID,MPI_COMM_WORLD,&nFpts_in);

  // Sloppy, but necessary to breakup communication from computation more efficiently
  isMPI = 1;
#endif
}

void mpiFace::finishRightSetup(void)
{
#ifndef _NO_MPI
  MPI_Wait(&nFpts_out,MPI_STATUSES_IGNORE);
  MPI_Wait(&nFpts_in,MPI_STATUSES_IGNORE);

  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up
   * NOTE THAT THIS IS DIFFERENT THAN IN intFaces --- */
  fptStartR = nFptsR;
  fptEndR = 0;

  UR.setup(nFptsR,nFields);
  bufUR.setup(nFptsR,nFields);
#endif
}

void mpiFace::communicate(void)
{
  getLeftState();
#ifndef _NO_MPI
  /* Send/Get data to/from right element [order reversed to match left ele] */
  MPI_Isend(UL.getData(),UL.getSize(),MPI_DOUBLE,procR,IDR,MPI_COMM_WORLD,&UL_out);
  MPI_Irecv(bufUR.getData(),UR.getSize(),MPI_DOUBLE,procR,ID,MPI_COMM_WORLD,&UR_in);
#endif
}

void mpiFace::getRightState(void)
{
#ifndef _NO_MPI
  // Make sure the communication is complete & transfer from buffer
  MPI_Wait(&UR_in,&status);
  if (params->viscous) MPI_Wait(&gradUR_in,&status);

  // Copy UR from the buffer to the proper matrix
  int fpt = 0;
  for (int i=fptStartR-1; i>=fptEndR; i--) {
    for (int j=0; j<nFields; j++) {
      UR(fpt,j) = bufUR(i,j);
    }

    fpt++;
  }
#endif
}

void mpiFace::setRightState(void)
{
#ifndef _NO_MPI
  // Right state handled by counterpart across boundary - do nothing.
#endif
}
