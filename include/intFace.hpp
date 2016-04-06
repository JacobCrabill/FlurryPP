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
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
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

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  //! Get data from the right element
  void getRightState(void);

  //! For viscous cases, get the solution gradient from the right element
  void getRightGradient(void);

  //! Put the calculated interface flux into the right element's memory
  void setRightStateFlux(void);

  //! Put the common solution into the right element's memory (viscous cases)
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Do nothing [not an inlet/outlet boundary]
  vector<double> computeMassFlux(void);

private:
  int faceID_R;              //! Right element's face ID
  int relRot;              //! Relative rotation of right element's face (for 3D)
  int fptStartR, fptEndR;
  vector<int> fptR;        //! Indices of flux points in right element

  bool isNew_R = true; //! Flag for initialization (esp. due to unblanking)

  /* --- Storage for all solution/geometry data at flux points [right state] --- */
  vector<matrix<double>> FR;   //! Flux array [nFpts, nDims, nFields]
  vector<double*> FnR;    //! Common normal flux for right ele [in ele's memory]
  Array<double*,2> dUcR;    //! Common solution for left ele (in ele's memory)  [nFpts, nFields]
  matrix<double> normR;   //! Unit outward normal at flux points  [nFpts, nDims]
  vector<double> dAR;     //! Local face-area equivalent at flux points
};
