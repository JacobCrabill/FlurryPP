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

#include <array>
#include <vector>

#include "global.hpp"

#include "geo.hpp"
#include "matrix.hpp"

struct tetra
{
  array<point,4> nodes;    //! Positions of nodes in tet
  //vector<point> qpts;      //! Locations of quadrature points in tet
  //vector<double> weights;  //! Weights associated with each quadrature point
  int donorID;             //! Donor-grid cell ID to which tetra belongs
};

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

  vector<point> qpts;       //! Locations of quadrature points in reference tetrahedron
  vector<double> weights;   //! Quadrature weights
  vector<double> shapeQpts; //! Values of tetrahedron shape basis at quadrature points

  /* !!! -------------------------------------------- !!! */
  /* use only one of these methods: the 'tetra' method or the Array/'big node list' version */

  vector<tetra> tets;  //! Tetrahedrons comprising the supermesh

  /* !!! -------------------------------------------- !!! */

  Array<point,2> tetNodes; //! Physical nodal positions of each tet in supermesh
  vector<int> tetDonorID;  //! Donor-grid cell ID for each tet

  /* !!! -------------------------------------------- !!! */

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
  vector<tetra> splitHexIntoTet(vector<point> &hexNodes);

  //! Use the given face and outward normal to clip the given tet and return the new set of tets
  vector<tetra> clipTet(tetra &tet, vector<point> &clipFace, Vec3 &norm);

};
