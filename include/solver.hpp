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

#include <map>
#include <set>

#include "global.hpp"
#include "ele.hpp"
#include "face.hpp"
#include "operators.hpp"
#include "input.hpp"

class solver
{
friend class geo; // Probably only needed if I make eles, opers private?

public:
  /* === Member Variables === */
  //! Map from eType to order to element- & order-specific operator
  map<int, map<int,oper> > opers;

  //! Vector of all eles handled by this solver
  vector<ele> eles;

  //! Vector of all faces handled by this solver
  vector<ele> faces;

  /* === Setup Functions === */
  solver();

  void initialize(input *params);

  void setupOperators();

  ~solver();

  /* === Functions Related to Basic FR Process === */

  //! Perform one full step of computation
  void calcResidual(void);

  //! Extrapolate the solution to the flux points
  void extrapolateU(void);

  //! Calculate the inviscid flux at the solution points
  void calcInviscidFlux_spts(void);

  //! Calculate the inviscid interface flux at all element faces
  void calcInviscidFlux_faces(void);

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

  //! Calculate the gradient of the flux at the solution points
  void calcGradF_spts(void);

  //! Calculate the divergence of the flux at the solution points
  void calcDivF_spts(void);

  //! Apply the correction function & add to the divergence of the flux
  void correctFlux(void);

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
  map<int,set<int>> polyOrders;

  //! Lists of cells to apply various adaptation methods to
  vector<int> r_adapt_cells, h_adapt_cells, p_adapt_cells;
}
