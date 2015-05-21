/*!
 * \file face.hpp
 * \brief Header file for the face class
 *
 * Class to handle calculation of interfce fluxes between elements
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
#pragma once

#include <vector>

#ifndef _NO_MPI
#include "mpi.h"
#endif

#include "matrix.hpp"
#include "face.hpp"

class face;

class mpiFace : public face
{
public:

  void setupRightState(void);

  /*! Wait to receive # of flux points from right side of proccessor */
  void finishRightSetup(void);

  void getRightState(void);

  void setRightState(void);

  /*! Perform the MPI communication across the processor boundary */
  void communicate(void);

  int procL;               //! Processor ID for left
  int procR;               //! Processor ID for right
  int IDR;                 //! Local face ID of face on right processor

private:
  int locF_R;              //! Right element's local face ID
  int fptStartR, fptEndR;

  /* --- Storage for all solution/geometry data at flux points [right state] --- */
  //vector<matrix<double>> FR;   //! Flux array [nFpts, nDims, nFields]
  matrix<double> normR;   //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacR; //! Determinant of transformation Jacobian

  matrix<double> bufUR;      //! Incoming buffer for receving UR
  matrix<double> bufGradUR;  //! Incoming buffer for receving gradUR

#ifndef _NO_MPI
  MPI_Request UL_out;
  MPI_Request UR_in;
  MPI_Request gradUL_out;
  MPI_Request gradUR_in;
  MPI_Request nFpts_out;
  MPI_Request nFpts_in;

  MPI_Status status;

  int* fptsBuffOut;
  int* fptsBuffIn;
#endif
};
