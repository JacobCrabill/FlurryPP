/*!
 * \file multigrid.cpp
 * \brief multigrid class definition
 *
 * The multigrid class stores a collection of solvers and grids, and handles
 * the communication between them to run cycles of h- or p-multigrid
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 */

#include "../include/multigrid.hpp"

#include "matrix.hpp"

multigrid::multigrid()
{

}

multigrid::~multigrid()
{

}

void multigrid::initialize(geo *Geo0, solver *Solver0, input* params)
{
  this->params = params;

  Solvers.resize(1);
  Geos.resize(1);

  Solvers[0] = Solver0;
  Geos[0] = Geo0;

  nGrids = 1;
}

void multigrid::setupNextFineLevel(int Lvl0)
{
  matrix<int> c2v = Geos[Lvl0]->c2v;
  vector<int> cell(4); //! currently assuming linear elements
  vector<int> icTmp(4);
  vector<point> xv = Geos[Lvl0]->xv;
  point pt;
  int nVerts = Geos[Lvl0]->nVerts;
  int nEles = Geos[Lvel0]->nEles;

  // match coarse-grid elements to fine-grid elements
  matrix<int> icC2F(nEles,4);
  vector<int> icF2C;

  for (int ic=0; ic < Geos[Lvl0]->nEles; ic++) {

    /* -- Split quad into tris -- */

    // Get midpoint of quad & add to xv
    pt.x = 0;  pt.y = 0;
    for (int i=0; i<4; i++) {
      pt.x += xv[c2v(ic,i)].x/4.;
      pt.y += xv[c2v(ic,i)].y/4.;
    }
    Geos[Lvl0]->xv.push_back(pt);

    // Insert new tris
    cell[0] = Geos[Lvl0]->c2v(ic,1);
    cell[1] = Geos[Lvl0]->c2v(ic,2);
    cell[2] = nVerts;
    cell[3] = nVerts;

    c2v.insertRow(cell);

    cell[0] = Geos[Lvl0]->c2v(ic,2);
    cell[1] = Geos[Lvl0]->c2v(ic,3);
    cell[2] = nVerts;
    cell[3] = nVerts;

    c2v.insertRow(cell);

    cell[0] = Geos[Lvl0]->c2v(ic,3);
    cell[1] = Geos[Lvl0]->c2v(ic,1);
    cell[2] = nVerts;
    cell[3] = nVerts;

    c2v.insertRow(cell);

    // Original quad->tri: leave first 2 points; make last 2 new center pt
    // (Duplicate center point for collapsed-edge tri)
    c2v(ic,2) = nVerts;
    c2v(ic,3) = nVerts;

    nVerts++;

    // Coarse to Fine, Fine to Coarse mapping
    icC2F(ic,0) = ic;    icF2C.push_back(ic);
    icC2F(ic,1) = nEles; icF2C.push_back(ic); nEles++;
    icC2F(ic,2) = nEles; icF2C.push_back(ic); nEles++;
    icC2F(ic,3) = nEles; icF2C.push_back(ic); nEles++;
  }

  /* Push back mappings into vector storage */
  mg_icC2F.push_back(icC2F);
  mg_icF2C.push_back(icF2C);

  mg_nEles.push_back(nEles);
  mg_nVerts.push_back(nVerts);

  /* Setup the new geo and solver objects */
  geo *Geo = new geo();
  Geo->c2v = c2v;
  Geo->xv = xv;
  Geo->bndPts = Geos[Lvl0]->bndPts;

  Geo->processConnectivity();
  Geo->processPeriodicBoundaries();

  solver *Solver = new solver();
  Solver->setup(params,Geo);

  Geos.push_back(Geo);
  Solvers.push_back(Solver);

  setupRestriction(Lvl0,nGrids);
  setupProlongation(Lvl0,nGrids);

  nGrids++;
}

void multigrid::restrictResidual(int coarse, int fine)
{

}

void multigrid::prolongResidual(int coarse, int fine)
{

}

void multigrid::setupRestriction(int coarse, int fine)
{
  for (int i=0; i<Solvers[coarsse]->eles.size(); i++) {

  }
}

void multigrid::setupProlongation(int coarse, int fine)
{
  for (int i=0; i<Solvers[fine]->eles.size(); i++) {
    vector<point> pts = Solvers[fine]->eles[i].getPosSpts();
    // get inverse mapping within Solvers[coarse]->eles[ic_F2C[i]]
    // Setup interpolation operator for this eles within coarse ele
  }
}
