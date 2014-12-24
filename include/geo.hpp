/*!
 * \file geo.hpp
 * \brief Header file for geometry class
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

#include "global.hpp"
#include "input.hpp"
#include "solver.hpp"

class geo
{
public:
  geo();

  /* === Primay setup routines === */

  //! Setup the geomery using input parameters
  setup(input* _params);

  //! Take the basic connectivity data and generate the rest
  processConnectivity();

  //! Create the elements and faces needed for the simulation
  setupElesFaces(solver* Solver);

  /* === Helper Routines === */

  //! Read essential connectivity from a Gmsh mesh file
  readGmsh();

  //! Create a simple Cartesian mesh from input parameters
  createMesh();

  createQuadMesh(); // or lump into createMesh()
  createTriMesh();

  ~geo();

private:

  input *params;
  int nDims, nFields;
  int nEles, nVerts, nFaces;

  // Basic [essential] Connectivity Data
  matrix<int> c2v;
  vector<point> xv;

  // Additional Connectivity Data
  matrix<int> c2e, e2c, e2v, v2e, v2v, v2nv, v2c, v2nc, c2nv, c2ne;
  vector<int> bndPts, bndFaces;
};
