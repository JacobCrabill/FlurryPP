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

#include "global.hpp"
#include "ele.hpp"
#include "matrix.hpp"
#include "input.hpp"

class face
{
public:
  face();

  /*! Setup access to the left & right elements' data */
  void setupFace(ele *eL, ele *eR, int locF_L, int locF_R, int gID);

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

  /* --- Storage for all solution/geometry data at flux points --- */
  //vector<vector<double>*> UL;
  //vector<vector<double>*> UL;
  vector<double*> UL, UR;         //! Discontinuous solution at left and right
  vector<double*> disFnL, disFnR; //! Discontinuous normal flux at left and right
  vector<matrix<double>*> gradUL;
  vector<matrix<double>*> gradUR;
  //vector<matrix<double>*> FL;
  //vector<matrix<double>*> FR;
  vector<matrix<double*>> FL; // will this even work...? if not, need vec<vec<vec<dbl*>>> ?
  vector<matrix<double*>> FR;
  //vector<matrix<double*>> FL, FR;
  //vector<vector<vector<double*>>> FL;
  //vector<vector<vector<double*>>> FR;
  //vector<vector<double>*> dFnL;   //! Common minus discontinuous normal flux for left ele
  //vector<vector<double>*> dFnR;   //! Common minus discontinuous normal flux for right ele
  vector<double*> dFnL;   //! Common minus discontinuous normal flux for left ele
  vector<double*> dFnR;   //! Common minus discontinuous normal flux for right ele
  matrix<double> Fn;     // Can't use ptr, b/c 2 eles - they need to point to this instead
  matrix<double> normL, normR; //! Unit outward normal at flux points
  vector<double> dAL, dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacL;
  vector<double> detJacR;

  matrix<double> tempFL, tempFR;
  vector<double> tempUL, tempUR;

  // Probably only needed for debugging... remove this later
  vector<point> posFpts;
};
