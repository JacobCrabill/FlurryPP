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

#include "global.hpp"

class ele;

#include "matrix.hpp"
#include "input.hpp"

class face
{
public:

  /*! Assign basic parameters to boundary */
  void initialize(ele *eL, ele *eR, int locF_L, int rightParam, int gID, input* params);

  /*! Setup access to the left elements' data */
  void setupFace(void);

  /*! Setup access to the right elements' data (if it exists) */
  virtual void setupRightState(void) =0;

  /*! Get the values of the solution to the left of the face */
  void getLeftState(void);

  /*! Get the values of the solution to the right of the face */
  virtual void getRightState(void) =0;

  /*! For all internal faces, put the normal flux into the right ele
   *  Either put directly into ele's memory, or send across MPI boundary */
  virtual void setRightState(void) =0;

  /*! Calculate the common inviscid flux on the face */
  void calcInviscidFlux(void);

  /*! Calculate the common viscous flux on the face */
  void calcViscousFlux(void);

  int ID; //! Global ID of face

  input *params; //! Input parameters for simulation

  double maxDU;
//protected:
  int nFptsL, nFptsR;
  int nDims, nFields;
  int locF_L;
  int fptStartL, fptEndL;
  int rightParam;

  ele* eL;
  ele* eR;

  /* --- Storage for all solution/geometry data at flux points [left state] --- */
  matrix<double> UL;      //! Discontinuous solution at left, right eles [nFpts, nFields]
  matrix<double> UR;      //! Discontinuous solution at left, right eles [nFpts, nFields]
  vector<matrix<double>> gradUL; //! Solution gradient at left side
  vector<matrix<double>> gradUR; //! Solution gradient at right side
  vector<matrix<double>> FL; //! Flux matrix at each flux point [nFpts, nDims, nFields]
  vector<double*> FnL;    //! Common normal flux for left ele (in ele's memory)  [nDims, nFpts, nFields]
  matrix<double> Fn;      //! Common numerical flux at interface  [nFpts, nFields]
  matrix<double> normL;   //! Unit outward normal at flux points
  vector<double> dAL;     //! Local face-area equivalent (aka edge Jacobian) at flux points
  vector<double> detJacL; //! Determinant of transformation Jacobian at flux points

  //! Temporary vectors for calculating common flux
  matrix<double> tempFL, tempFR;
  vector<double> tempUL;

  // Probably only needed for debugging... remove this later
  vector<point> posFpts;

  int isMPI;  //! Flag for MPI faces to separate communication from flux calculation
};
