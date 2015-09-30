/*!
 * \file boundFace.hpp
 * \brief Header file for the boundFace class
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

  /*! Assign boundary-condition type */
  void setupRightState(void);

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  /*! Apply invisic boundary conditions to the solution */
  void applyBCs(void);

  /*! Apply viscous boundary conditions to the solution */
  void applyViscousBCs(void);

  /*! Call applyBCs */
  void getRightState(void);

  /*! Call applyBCs */
  void getRightGradient(void);

  /*! No right element at a boundary - do nothing. */
  void setRightStateFlux(void);

  /*! No right element at a boundary - do nothing. */
  void setRightStateSolution(void);

  /*! For wall boundary conditions, compute the force on the wall */
  vector<double> computeWallForce(void);

private:
  int bcType;  //! Boundary condition to apply to this face

  matrix<double> deltaU;
  matrix<double> deltaUdot;
  matrix<double> deltaUint;
};
