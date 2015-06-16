/*!
 * \file mpiFace.cpp
 * \brief Class to handle interface flux calculations between processors
 *
 * Handles MPI communication to its counterpart across an MPI boundary
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

  /* Send/Get # of flux points to/from right element */
  MPI_Irecv(&nFptsR,1,MPI_INT,procR,ID,MPI_COMM_WORLD,&nFpts_in);
  MPI_Isend(&nFptsL,1,MPI_INT,procR,IDR,MPI_COMM_WORLD,&nFpts_out);

  // Sloppy, but necessary to breakup communication from computation more efficiently
  isMPI = 1;
#endif
}

void mpiFace::finishRightSetup(void)
{
#ifndef _NO_MPI
  MPI_Wait(&nFpts_in,MPI_STATUSES_IGNORE);
  MPI_Wait(&nFpts_out,MPI_STATUSES_IGNORE);  

  /* --- Will have to introduce 'mortar' elements in the future [for p-adaptation],
   * but for now just force all faces to have same # of flux points [order] --- */

  if (nFptsL != nFptsR)
    FatalError("Mortar elements not yet implemented - must have nFptsL==nFptsR");

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up ---
   * NOTE THAT THIS IS DIFFERENT THAN IN intFaces */
  fptStartR = nFptsR;
  fptEndR = 0;

  UR.setup(nFptsR,nFields);
  bufUR.setup(nFptsR,nFields);
  bufGradUR.setup(nFptsR,nDims*nFields); // !! TEMP HACK !!
  bufGradUL.setup(nFptsR,nDims*nFields); // !! TEMP HACK !!
#endif
}

void mpiFace::communicate(void)
{
  getLeftState();

#ifndef _NO_MPI
  /* Send/Get data to/from right element [order reversed to match left ele] */

  // The send/receive pairs are tagged by the processor-local face ID of the
  // face on the receiving end of the call
  MPI_Irecv(bufUR.getData(),UR.getSize(),MPI_DOUBLE,procR,ID,MPI_COMM_WORLD,&UR_in);
  MPI_Isend(UL.getData(),UL.getSize(),MPI_DOUBLE,procR,IDR,MPI_COMM_WORLD,&UL_out);

  if (params->viscous) {
    // !!! TEMP HACK !!! Just until I update Matrix class to 3D+
    for (int i=0; i<nFptsL; i++)
      for (int j=0; j<nDims; j++)
        for (int k=0; k<nFields; k++)
          bufGradUL(i,j+k*nDims) = gradUL[i](j,k);

    MPI_Irecv(bufGradUR.getData(),bufGradUR.getSize(),MPI_DOUBLE,procR,ID,MPI_COMM_WORLD,&gradUR_in);
    MPI_Isend(bufGradUL.getData(),bufGradUL.getSize(),MPI_DOUBLE,procR,IDR,MPI_COMM_WORLD,&gradUL_out);
  }
#endif
}

void mpiFace::getRightState(void)
{
#ifndef _NO_MPI
  // Make sure the communication is complete & transfer from buffer
  MPI_Wait(&UL_out,&status);
  MPI_Wait(&UR_in,&status);
  if (params->viscous) {
    MPI_Wait(&gradUL_out,&status);
    MPI_Wait(&gradUR_in,&status);
  }

  // Copy UR from the buffer to the proper matrix [note that the order of the
  // fpts is reversed between the two faces]
  int fpt = 0;
  for (int i=fptStartR-1; i>=fptEndR; i--) {
    for (int j=0; j<nFields; j++)
      UR(fpt,j) = bufUR(i,j);

    if (params->viscous) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUR[fpt](dim,j) = bufUR(i,dim+j*nDims);
    }

    fpt++;
  }
#endif
}

void mpiFace::setRightStateFlux(void)
{
#ifndef _NO_MPI
  // Right state handled by counterpart across boundary - do nothing.
#endif
}

void mpiFace::setRightStateSolution(void)
{
#ifndef _NO_MPI
  // Right state handled by counterpart across boundary - do nothing.
#endif
}
