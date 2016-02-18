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

  /*! For inlet/outlet boundary conditions, compute the force on the wall */
  vector<double> computeMassFlux(void);

private:
  int bcType;  //! Boundary condition to apply to this face

  matrix<double> deltaU;
  matrix<double> deltaUdot;
  matrix<double> deltaUint;
};
