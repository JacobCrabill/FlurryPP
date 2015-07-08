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
  IDR = myInfo.IDR;
  relRot = myInfo.relRot;
  procL = myInfo.procL;
  procR = myInfo.procR;
  myComm = myInfo.gridComm;

  /* Send/Get # of flux points to/from right element */
  MPI_Irecv(&nFptsR,1,MPI_INT,procR,ID,*myComm,&nFpts_in);
  MPI_Isend(&nFptsL,1,MPI_INT,procR,IDR,*myComm,&nFpts_out);

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

  /* --- Setup the L/R flux-point matching ---
   * NOTE THAT THIS IS DIFFERENT THAN IN intFaces */
  fptR.resize(nFptsL);
  if (nDims == 2) {
    // For 1D faces [line segments] only - find first/last ID of fpts;
    // right faces's points are simply reversed
    fptStartR = nFptsR;
    fptEndR = 0;

    int fpt = 0;
    for (int i=fptStartR-1; i>=fptEndR; i--) {
      fptR[fpt] = i;
      fpt++;
    }
  }
  else if (nDims == 3) {
    // Only for quad tensor-product faces: Rotate the face to the correct relative orientation
    int order = sqrt(nFptsL)-1;
    for (int i=0; i<nFptsL; i++) {
      int ifpt = i%(order+1);
      int jfpt = floor(i/(order+1));
      switch (relRot) {
        case 0:
          fptR[i] = ifpt*(order+1) + jfpt;
          break;
        case 1:
          fptR[i] = order-ifpt + jfpt*(order+1);
          break;
        case 2:
          fptR[i] = nFptsL-1 - (ifpt*(order+1) + jfpt);
          break;
        case 3:
          fptR[i] = nFptsL - (order+1)*(jfpt+1) + ifpt;
          break;
      }
    }
  }

  UR.setup(nFptsR,nFields);
  bufUR.setup(nFptsR,nFields);
  bufGradUR.setup(nFptsR,nDims*nFields); // !! TEMP HACK !!  need 3D matrix/array
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
  MPI_Irecv(bufUR.getData(),UR.getSize(),MPI_DOUBLE,procR,ID,*myComm,&UR_in);
  MPI_Isend(UL.getData(),UL.getSize(),MPI_DOUBLE,procR,IDR,*myComm,&UL_out);

  if (params->viscous) {
    // !!! TEMP HACK !!! Just until I update Matrix class to 3D+
    for (int i=0; i<nFptsL; i++)
      for (int j=0; j<nDims; j++)
        for (int k=0; k<nFields; k++)
          bufGradUL(i,j+k*nDims) = gradUL[i](j,k);

    MPI_Irecv(bufGradUR.getData(),bufGradUR.getSize(),MPI_DOUBLE,procR,ID,*myComm,&gradUR_in);
    MPI_Isend(bufGradUL.getData(),bufGradUL.getSize(),MPI_DOUBLE,procR,IDR,*myComm,&gradUL_out);
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
  for (int i=0; i<nFptsL; i++) {
    for (int j=0; j<nFields; j++)
      UR(fpt,j) = bufUR(fptR[i],j);

    if (params->viscous) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUR[fpt](dim,j) = bufGradUR(fptR[i],dim+j*nDims);
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
