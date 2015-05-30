/*!
 * \file multigrid.hpp
 * \brief Header file for multigrid class
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

#include "solver.hpp"
#include "geo.hpp"
#include "global.hpp"
#include "input.hpp"

class multigrid
{
public:

  /* --- Member Variables --- */

  vector<solver*> Solvers;
  vector<geo*> Geos;

  vector<matrix<int>> mg_icC2F; //! Mappings from coarse grid cells to fine grid cells
  vector<vector<int>> mg_icF2C; //! Mapping from fine grid cells to coarse grid cells
  vector<int> mg_nEles;
  vector<int> mg_nVerts;

  input* params;

  int nGrids;

  /* --- Member Functions --- */

  multigrid();
  ~multigrid();

  //! Initialize the multigrid
  void initialize(solver* Solver0, geo* Geo0, input* params);

  //! Initialize the multigrid
  void setLevel0(geo* Geo0, solver* Solver0);

  //! Use an existing grid level to setup another
  void setupNextLevel(int level_0);

  //! Setup the next-coarser grid level from the given level
  void setupNextCoarseLevel(int level_0);

  //! Setup the next-finer grid level from the given level
  void setupNextFineLevel(int Lvl0);

  //! Restrict the current solution residual from one grid to another
  void restrictResidual(int coarse, int fine);

  //! Restrict the solution update from one grid to another
  void prolongResidual(int coarse, int fine);

  //! Setup the prolongation operator between the two grids
  void setupProlongation(int coarse, int fine);

  //! Setup the restriction operator between the two grids
  void setupRestriction(int coarse, int fine);
};
