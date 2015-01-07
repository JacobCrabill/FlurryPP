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

#include "global.hpp"
#include "ele.hpp"
#include "matrix.hpp"
#include "input.hpp"

class face
{
public:
  face();

  /*! Setup access to the left & right elements' data */
  setupFace(ele *eL, ele *eR, int locF_L, int locF_R, int gID);

  /*! Calculate the common inviscid flux on the face */
  calcInviscidFlux();

  /*! Calculate the common viscous flux on the face */
  calcViscousFlux();

  ~face();

  int ID; //! Global ID of face

private:
  int nFptsL, nFptsR;
  int nDims, nFields;

  input *params;

  /* --- Storage for all solution/geometry data at flux points --- */
  vector<vector<double>*> UL;
  vector<vector<double>*> UR;
  vector<matrix<double>*> gradUL;
  vector<matrix<double>*> gradUR;
  vector<matrix<double>*> FL;
  vector<matrix<double>*> FR;
  vector<vector<double>*> Fn;
  vector<vector<double>*> deltaF;
  vector<array<double,3>*> normL;
  vector<array<double,3>*> normR;

  matrix<double> tempFL, tempFR;
  vector<double> tempUL, tempUR;
};
