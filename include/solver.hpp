/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include <map>
#include <set>
#include <vector>

class oper;

#include "global.hpp"

#include "bound.hpp"
#include "ele.hpp"
#include "face.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "operators.hpp"

class solver
{
friend class geo; // Probably only needed if I make eles, opers private?

public:
  /* === Member Variables === */
  //! Pointer to geometry object for mesh-related operations
  geo *Geo;

  //! Map from eType to order to element- & order-specific operator
  map<int, map<int,oper> > opers;

  //! Vector of all eles handled by this solver
  vector<ele> eles;

  //! Vector of all interior faces handled by this solver
  vector<face> faces;

  //! Vector of all boundary faces handled by this solver
  vector<bound> bounds;

  /* === Setup Functions === */
  solver();

  void setup(input *params, geo *Geo);

  //! Setup the FR operators for all ele types and polynomial orders which will be used in computation
  void setupOperators();

  /* === Functions Related to Basic FR Process === */

  //! Apply the initial condition to all elements
  void initializeSolution();

  void update(void);

  //! Perform one full step of computation
  void calcResidual(int step);

  //! Advance solution in time
  void timeStepA(int step);

  void timeStepB(int step);

  void copyUspts_U0(void);
  void copyU0_Uspts(void);

  //! Extrapolate the solution to the flux points
  void extrapolateU(void);

  //! Extrapolate the solution to the mesh (corner) points
  void extrapolateUMpts(void);

  //! Calculate the inviscid flux at the solution points
  void calcInviscidFlux_spts(void);

  //! Calculate the inviscid interface flux at all element faces
  void calcInviscidFlux_faces(void);

  //! Calculate the inviscid interface flux at all boundary faces
  void calcInviscidFlux_bounds(void);

  //! Calculate the gradient of the solution at the solution points
  void calcGradU_spts(void);

  //! For viscous calculations, apply the correction procedure to the solution
  void correctU(void);

  /*! For viscous calculations, extrapolate the corrected gradient of the solution
   *  from the solution points to the flux points */
  void extrapolateGradU(void);

  //! Calculate the viscous flux at the solution points
  void calcViscousFlux_spts(void);

  //! Calculate the viscous interface flux at all element faces
  void calcViscousFlux_faces(void);

  //! Calculate the viscous interface flux at all boundary faces
  void calcViscousFlux_bounds(void);

  //! Calculate the gradient of the flux at the solution points
  void calcGradF_spts(void);

  //! Use full space-time chain-rule to transform gradients to parent domain
  void transformGradF_spts(void);

  //! Calculate the divergence of the flux at the solution points
  void calcDivF_spts(int step);

  //! Extrapolate total flux to flux points & dot with normal
  void extrapolateNormalFlux(void);

  //! Apply the correction function & add to the divergence of the flux
  void correctDivFlux(int step);

  // **All of the following functions are just food for thought at the moment**

  /* === Functions for Shock Capturing & Filtering=== */
  void shockCapturing();

  /* === Functions Related to Adaptation === */
  void get_r_adapt_cells();

  void get_p_adapt_cells();

  void get_h_adapt_cells();

  void setup_r_adaption();

  void setup_h_adaptation();

  void setup_p_adaptation();

  void add_ele(int eType, int order);

  /* === Functions Related to Overset Grids === */

private:
  //! Pointer to the parameters object for the current solution
  input *params;

  //! Set of all element types present in current simulation
  set<int> eTypes;

  //! Set of all polynomial orders present for each element type
  map<int,set<int> > polyOrders;

  //! Lists of cells to apply various adaptation methods to
  vector<int> r_adapt_cells, h_adapt_cells, p_adapt_cells;

  int nRKSteps;

  vector<double> RKa, RKb;
};
