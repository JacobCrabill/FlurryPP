/*!
 * \file intFace.hpp
 * \brief Header file for the intFace class
 *
 * Class to handle calculation of interfce fluxes between all elements
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

class face;

#include "face.hpp"

class intFace : public face
{
public:
  //! Setup arrays to handle getting/setting data at right element
  void setupRightState(void);

  //! Get data from the right element
  void getRightState(void);

  //! Put the calculated interface flux into the right element's memory
  void setRightStateFlux(void);

  //! Put the common solution into the right element's memory (viscous cases)
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

private:
  int locF_R;              //! Right element's local face ID
  int relRot;              //! Relative rotation of right element's face (for 3D)
  int fptStartR, fptEndR;
  vector<int> fptR;        //! Indices of flux points in right element

  /* --- Storage for all solution/geometry data at flux points [right state] --- */
  vector<matrix<double>> FR;   //! Flux array [nFpts, nDims, nFields]
  vector<double*> FnR;    //! Common normal flux for right ele [in ele's memory]
  vector<double*> UcR;    //! Common solution for left ele (in ele's memory)  [nFpts, nFields]
  matrix<double> normR;   //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;     //! Local face-area equivalent at flux points
  vector<double> detJacR; //! Determinant of transformation Jacobian
};
