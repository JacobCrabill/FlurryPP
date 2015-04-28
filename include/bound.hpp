/*!
 * \file bound.hpp
 * \brief Header file for the bound class
 *
 * Class to handle enforcement of boundary conditions at boundary faces
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

#include "global.hpp"
#include "input.hpp"
#include "ele.hpp"
#include "flux.hpp"

class bound
{
public:

  /*! Setup access to the left & right elements' data */
  void setupBound(ele *eL, int locF_L, int bcType, int gID);

  /*! Calculate the inviscid flux from the boundary condition */
  void calcInviscidFlux(void);

  /*! Calculate the viscous flux from the boundary condition */
  void calcViscousFlux(void);

  int ID; //! Global ID of face

  input *params; //! Input parameters for simulation

  void applyBCs(const double *uL, double* uR, const double* norm);

private:
  int nFptsL;
  int nDims, nFields;
  int bcType;
  int locF_L;

  /* --- Storage for all solution/geometry data at flux points --- */
  vector<double*> UL;       //! Discontinuous solution at element boundary
  matrix<double> UR;        //! Boundary condition from "ghost right state"
  vector<double*> disFnL;   //! Discontinuous normal flux at element boundary
  vector<matrix<double>*> gradUL;
  vector<matrix<double*>> FL;
  vector<double*> dFnL;     //! Common minus discontinuous normal flux for ele
  matrix<double> Fn;        //! Common normal flux on boundary
  vector<double*> deltaF;
  matrix<double> normL;
  vector<double> dAL;       //! Local face-area equivalent at flux points
  vector<double> detJacL;

  matrix<double> tempFL, tempFR;
  vector<double> tempUL;

  vector<vector<double>> UC;
};
