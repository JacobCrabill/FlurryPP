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

#include "ele.hpp"
#include "matrix.hpp"
#include "input.hpp"

class face
{
public:
  face();

  /*! Setup access to the left & right elements' data */
  void setupFace(ele *eL, ele *eR, int locF_L, int locF_R, int gID, input* params);

  /*! Calculate the common inviscid flux on the face */
  void calcInviscidFlux(void);

  /*! Calculate the common viscous flux on the face */
  void calcViscousFlux(void);

  int ID; //! Global ID of face

  input *params; //! Input parameters for simulation

private:
  int nFptsL, nFptsR;
  int nDims, nFields;
  int locF_L, locF_R;

  ele* eL;
  ele* eR;

  /* --- Storage for all solution/geometry data at flux points --- */
  vector<double*> UL, UR;         //! Discontinuous solution at left and right
  vector<matrix<double>*> gradUL;
  vector<matrix<double>*> gradUR;
  vector<matrix<double*>> FL; // will this even work...? if not, need vec<vec<vec<dbl*>>> ?
  vector<matrix<double*>> FR;
  vector<double*> FnL;   //! Common minus discontinuous normal flux for left ele
  vector<double*> FnR;   //! Common minus discontinuous normal flux for right ele
  matrix<double> Fn;     // Can't use ptr, b/c 2 eles - they need to point to this instead
  vector<double*> normL, normR; //! Unit outward normal at flux points
  vector<double*> dAL, dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacL;
  vector<double> detJacR;

  matrix<double> tempFL, tempFR;
  vector<double> tempUL, tempUR;

  // Probably only needed for debugging... remove this later
  vector<point> posFpts;
};
