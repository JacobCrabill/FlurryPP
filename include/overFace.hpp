/*!
 * \file overFace.hpp
 * \brief Header file for the overFace class
 *
 * Class to handle calculation of interfce fluxes on overset faces
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */
#pragma once

class face;
class overComm;

#include "face.hpp"

#include "overComm.hpp"

class overFace : public face
{
public:
  //! Overset communicator object where the face's data will be found
  overComm *OComm;

  //! Setup arrays to handle getting/setting data from Solver
  void setupRightState(void);

  //! Get the interpolated overset data from the Solver
  void getRightState(void);

  //! Do nothing [right state is non-existant]
  void setRightStateFlux(void);

  //! Do nothing [right state is non-existant]
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Return the physical position of the face's flux points
  vector<point> getPosFpts(void);

  int fptOffset;         //! Offset within Solver's mesh block-global interp point list
  vector<point> posFpts; //! Physical locations of left ele's flux points

private:

};
