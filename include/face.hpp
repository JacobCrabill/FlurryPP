/*!
 * \file face.hpp
 * \brief Header file for the abstract face parent class
 *
 * Abstract parent class to handle calculation of interfce fluxes at all faces
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <memory>
#include <vector>

#ifndef _NO_MPI
#include "mpi.h"
#endif

#include "global.hpp"

class ele;

#include "matrix.hpp"
#include "input.hpp"

//! Struct to pass into each face; values assigned as needed based on face type
struct faceInfo {
  int isMPI = 0;  //! Flag for whether face is an MPI-boundary face
  int isBnd = 0;  //! Flag for whether face is a boundary-condition face
  int IDR;

  //! Relative orientation between 3D faces
  int relRot;

  //! boundFace parameters
  int bcType;

  //! mpiFace parameters
  int procL;
  int procR;

#ifndef _NO_MPI
  MPI_Comm gridComm; //! MPI communicator for this grid
#endif
};

class face
{
public:

  /*! Assign basic parameters to boundary */
  void initialize(shared_ptr<ele> &eL, shared_ptr<ele> &eR, int gID, int locF_L, struct faceInfo myInfo, input* params);

  /*! Setup arrays and access to the left elements' data */
  void setupFace(void);

  /*! Setup pointer access to left elements' data */
  void getPointers(void);

  /*! Get pointer access to right element's data */
  virtual void getPointersRight(void) =0;

  /*! Setup access to the right elements' data (if it exists) */
  virtual void setupRightState(void) =0;

  /*! Get the values of the solution to the left of the face */
  void getLeftState(void);

  /*! For viscous cases, get solution gradient to the left of the face */
  void getLeftGradient(void);

  /*! Get the values of the solution to the right of the face */
  virtual void getRightState(void) =0;

  /*! For viscous cases, get the solution gradient to the right of the face */
  virtual void getRightGradient(void) =0;

  /*! For all internal faces, put the normal flux into the right ele
   *  Either put directly into ele's memory, or send across MPI boundary */
  virtual void setRightStateFlux(void) =0;

  /*! Viscous cases: For all internal faces, put the common solution into the right ele
   *  Either put directly into ele's memory, or send across MPI boundary */
  virtual void setRightStateSolution(void) =0;

  /*! Compute the force on any wall boundary conditions */
  virtual vector<double> computeWallForce(void) =0;

  /*! Calculate the common inviscid flux on the face */
  void calcInviscidFlux(void);

  /*! Calculate the common viscous flux on the face */
  void calcViscousFlux(void);

  /*! Calculate the common flux using the Rusanov method
   *  NOTE: Overridden for overFaces due to flux-interp modifications */
  virtual void rusanovFlux(void);

  /*! Calculate the common flux using the Roe method */
  void roeFlux(void);

  /*! Calculate the common flux using the Lax-Friedrichs method [scalar advection] */
  void laxFriedrichsFlux(void);

  /*! For boundary faces, use a central flux (no added dissipation) */
  void centralFluxBound(void);

  /*! Calculate a biased-average solution for LDG viscous flux */
  void ldgSolution(void);

  int ID; //! Global ID of face

  input *params; //! Input parameters for simulation

  int nFptsL, nFptsR;
  int nDims, nFields;
  int locF_L;
  int fptStartL, fptEndL;
  vector<int> rightParams;
  struct faceInfo myInfo;

protected:
  shared_ptr<ele> eL;
  shared_ptr<ele> eR;

  /* --- Storage for all solution/geometry data at flux points [left state] --- */
  matrix<double> UL;      //! Discontinuous solution at left ele [nFpts, nFields]
  matrix<double> UR;      //! Discontinuous solution at right ele [nFpts, nFields]
  matrix<double> UC;      //! Common solution at interface [nFpts, nFields]
  vector<matrix<double>> gradUL; //! Solution gradient at left side
  vector<matrix<double>> gradUR; //! Solution gradient at right side
  //Array<double,3> gradUR; //! Solution gradient at right side
  matrix<double> Vg;      //! Grid velocity at interface
  vector<matrix<double>> FL; //! Flux matrix at each flux point [nFpts, nDims, nFields]
  vector<double*> FnL;    //! Common normal flux for left ele (in ele's memory)  [nFpts, nFields]
  vector<double*> UcL;    //! Common solution for left ele (in ele's memory)  [nFpts, nFields]
  matrix<double> Fn;      //! Common numerical flux at interface  [nFpts, nFields]
  matrix<double> normL;   //! Unit outward normal at flux points
  vector<double> dAL;     //! Local face-area equivalent (aka edge Jacobian) at flux points
  vector<double> detJacL; //! Determinant of transformation Jacobian at flux points
  vector<double*> waveSp; //! Maximum numerical wave speed at flux point (in left ele's memory)

  //! Temporary vectors for calculating common flux
  matrix<double> tempFL, tempFR;
  vector<double> tempUL;

  int isMPI;  //! Flag for MPI faces to separate communication from flux calculation
  int isBnd;  //! Flag for boundary faces for use in LDG routines
};
