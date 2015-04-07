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

  void applyBCs(void);

private:
  int nFptsL;
  int nDims, nFields;
  int bcType;
  int locF_L;

  /* --- Storage for all solution/geometry data at flux points --- */
  vector<double*> UL;
  vector<matrix<double>*> gradUL;
  vector<matrix<double*>> FL;
  vector<double*> Fn;
  vector<double*> deltaF;
  matrix<double> normL;
  vector<double> detJacL;

  matrix<double> tempFL;
  vector<double> tempUL;

  vector<vector<double>> UC;
};
