/*!
 * \file superMesh.hpp
 * \brief Header file for superMesh class
 *
 * Creates a local supermesh for an element in one grid from
 * elements in another grid (See Farrell and Maddison, 2010)
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

#include <vector>

#include "global.hpp"

#include "geo.hpp"
#include "matrix.hpp"

class superMesh
{
public:

  /* ---- Default Functions ---- */

  superMesh();
  ~superMesh();

  /* ---- Member Variables ---- */

  geo* gridT;              //! Target grid which contains the target cell
  geo* gridD;              //! Donor grid in which to find the donor cells

  int targetCell;          //! Target cell ID for which to create local supermesh
  vector<int> donorCells;  //! Donor cell IDs from donor grid which overlap target cell

  int nTets;               //! Total number of tets comprising the supermesh
  int nQpts;               //! Total number of quadrature points in the whole supermesh
  int order;               //! Order of quadrature rule to use
  int nQpts_tet;           //! Number of quadrature points per tet (based on order)

  vector<point> locQpts;    //! Locations of quadrature points in reference tetrahedron
  vector<double> shapeQpts; //! Values of tetrahedron shape basis at quadrature points

  /* ---- Member Functions ---- */

  void setup(geo* _gridT, geo* gridD, int _targetCell, int _order);

  //! Using given grids and target cell, build the local supermesh
  void buildSuperMesh(void);

  //! Integrate a quantity over the superMesh (given at quadrature points of supermesh tets)
  double integrate(vector<double> &data);

  /*!
   * Returns the physical(?) positions of the quadrature points, along
   * with which donor-grid cell they lie within
   */
  void getQpts(vector<point> &qptPos, vector<int> &qptCell);

private:

  //! Subdivide the given hexahedron into 5 tetrahedrons
  Array<double,3> splitHexIntoTet(matrix<double> &hexNodes);

  //! Use the given face and outward normal to clip the given tet and return the new set of tets
  Array<double,3> clipTet(matrix<double> &tetNodes, matrix<double> &clipFace, Vec3 &norm);

};
