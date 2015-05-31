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

class face;

#include "face.hpp"

class boundFace : public face
{
public:

  void setupRightState(void);

  //! Apply boundary conditions to the solution
  void applyBCs(const double *uL, double* uR, const double* norm, int fpt);

  void getRightState(void);

  void setRightState(void);

//private:
  int bcType;  //! Boundary condition to apply to this face

  matrix<double> deltaU;
  matrix<double> deltaUdot;
  matrix<double> deltaUint;
};
