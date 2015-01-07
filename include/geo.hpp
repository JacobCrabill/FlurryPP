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
  void setup(input* params);

  //! Take the basic connectivity data and generate the rest
  void processConnectivity();

  //! Create the elements and faces needed for the simulation
  void setupElesFaces(solver* Solver);

  /* === Helper Routines === */

  //! Read essential connectivity from a Gmsh mesh file
  void readGmsh(string fileName);

  //! Create a simple Cartesian mesh from input parameters
  void createMesh();

  void createQuadMesh(); // or lump into createMesh()
  void createTriMesh();

  //! Get the reference-domain location of the solution points for the given element & polynomial order
  vector<point> getLocSpts(int eType, int order);

  //! Get the reference-domain location of the flux points for the given element & polynomial order
  vector<point> getLocFpts(int eType, int order);

  //! For tensor-product elements, get location of solution points along each direction
  vector<double> getLocSpts1D(int eType, int order);

  //! Get the 1D locations of flux points (for 2D elements and 3D tensor-product elements)
  vector<double> getLocFpts1D(int eType, int order);

  int nDims, nFields;
  int nEles, nVerts, nEdges, nFaces;

private:

  input *params;

  // Basic [essential] Connectivity Data
  matrix<int> c2v;
  vector<point> xv;

  // Additional Connectivity Data
  matrix<int> c2e, e2c, e2v, v2e, v2v, v2c;
  vector<int> v2nv, v2nc, c2nv, c2ne;
  vector<int> bndPts, bndFaces;
};
