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

  /*! Get pointer access to right ele's data */
  void getPointersRight(void);

  //! Get the interpolated overset solution data from the Solver
  void getRightState(void);

  //! Get the interpolated overset gradient data from the Solver
  void getRightGradient(void);

  //! Do nothing [right state is non-existant]
  void setRightStateFlux(void);

  //! Do nothing [right state is non-existant]
  void setRightStateSolution(void);

  //! Do nothing [not a wall boundary]
  vector<double> computeWallForce(void);

  //! Do nothing [not an inlet/outlet boundary]
  vector<double> computeMassFlux(void);

  //! Return the physical position of the face's flux points
  vector<point> getPosFpts(void);

  //! Return the outward unit normal at the face's flux points
  vector<point> getNormFpts();

  //! Override normal version when using flux-interp method
  void rusanovFlux(void);

  int fptOffset;         //! Offset within Solver's mesh block-global interp point list
  vector<point> posFpts; //! Physical locations of left ele's flux points

private:

};
