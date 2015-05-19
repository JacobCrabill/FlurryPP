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

#include "matrix.hpp"
#include "face.hpp"

class face;

class mpiFace : public face
{
public:

  void setupRightState(void);

  void getRightState(void);

  void setRightState(void);

private:
  int procL;               //! Processor ID for left
  int procR;               //! Processor ID for right
  int idR;                 //! Local face ID of face on right processor
  int locF_R;              //! Right element's local face ID
  int fptStartR, fptEndR;

  /* --- Storage for all solution/geometry data at flux points [right state] --- */
  vector<matrix<double>> gradUR;
  vector<matrix<double>> FR;   //! Flux array [nFpts, nDims, nFields]
  vector<double*> FnR;    //! Common normal flux for right ele [in ele's memory]
  matrix<double> normR;   //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacR; //! Determinant of transformation Jacobian
};
