/*!
 * \file solver_overset.cpp
 * \brief Overset-related methods for solver class
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

#include "solver.hpp"

#include <sstream>
#include <omp.h>

#include "input.hpp"
#include "geo.hpp"

/* ---- My New Overset Grid Functions ---- */

void solver::oversetInterp(void)
{
  OComm->exchangeOversetData(eles,opers);
}

void solver::setupOverset(void)
{
  if (gridRank == 0) cout << "Solver: Grid " << gridID << ": Setting up overset connectivity" << endl;

  OComm = make_shared<overComm>();

  OComm->setup(params,nGrids,gridID,gridRank,nprocPerGrid,Geo->gridIdList);

  OComm->matchOversetPoints(eles,overFaces,Geo->c2ac,Geo->eleMap,Geo->centroid,Geo->extents);
}

void solver::updateOversetConnectivity(void)
{
  Geo->updateOversetConnectivity();

  Geo->setupUnblankElesFaces(eles,faces,mpiFaces,overFaces);

  OComm->matchOversetPoints(eles,overFaces,Geo->c2ac,Geo->eleMap,Geo->centroid,Geo->extents);

  OComm->matchUnblankCells(eles,Geo->unblankCells,Geo->c2c,Geo->eleMap,params->order);
}

/* ---- Basic Tioga-Based Overset-Grid Functions ---- */

//void solver::setupOverset(void)
//{
//  if (gridRank == 0) cout << "Solver: Grid " << gridID << ": Setting up overset connectivity" << endl;
//  setupOversetData();
//  // Give TIOGA a pointer to this solver for access to callback functions
//  tg->setcallback(this);
//  Geo->updateOversetConnectivity();
//}

void solver::setupOversetData(void)
{
  // Allocate storage for global solution vector (for use with Tioga)
  // and initialize to 0
  int nSptsTotal = 0;
  for (uint i=0; i<eles.size(); i++) nSptsTotal += eles[i].getNSpts();
  U_spts.assign(nSptsTotal*params->nFields,0);
}

void solver::setGlobalSolutionArray(void)
{
  int ind = 0;
  for (uint i=0; i<eles.size(); i++) {
    eles[i].getUSpts(&U_spts[ind]);
    ind += eles[i].getNSpts() * params->nFields;
  }
}

void solver::updateElesSolutionArrays(void)
{
  int ind = 0;
  for (uint i=0; i<eles.size(); i++) {
    eles[i].setUSpts(&U_spts[ind]);
    ind += eles[i].getNSpts() * params->nFields;
  }
}

void solver::callDataUpdateTIOGA(void)
{
#ifndef _NO_MPI
  //cout << "Calling dataUpdate_highorder" << endl;
  tg->dataUpdate_highorder(params->nFields,U_spts.data(),0);
#endif
}

/* ---- Callback Functions for TIOGA Overset Grid Library ---- */

void solver::getNodesPerCell(int* cellID, int* nNodes)
{
  (*nNodes) = eles[*cellID].getNSpts();
}

void solver::getReceptorNodes(int* cellID, int* nNodes, double* posNodes)
{
  if (*cellID >= eles.size()) cout << "Invalid cellID!  cellID = " << *cellID << endl;
  eles[*cellID].getPosSpts(posNodes);
}

void solver::donorInclusionTest(int* cellID, double* xyz, int* passFlag, double* rst)
{
  // Determine if point is in cell: [x,y,z] should all lie between [-1,1]
  point refPt;
  (*passFlag) = eles[*cellID].getRefLocNelderMeade(point(xyz),refPt);

  rst[0] = refPt.x;
  rst[1] = refPt.y;
  rst[2] = refPt.z;
}

void solver::donorWeights(int* cellID, double* xyz, int* nWeights, int* iNode, double* weights, double* rst, int* fracSize)
{
  int ic = *cellID;

  // Get starting offset for global solution-point index
  int iStart = 0;
  for (int i=0; i<ic; i++)
    iStart += eles[i].getNSpts();

  // Put the global indices of ele's spts into inode array
  (*nWeights) = eles[ic].getNSpts();
  for (int i=0; i<(*nWeights); i++)
    iNode[i] = iStart + i;

  opers[eles[ic].eType][eles[ic].order].getInterpWeights(rst,weights);
}

void solver::convertToModal(int* cellID, int* nPtsIn, double* uIn, int* nPtsOut, int* iStart, double* uOut)
{
  // This is apparently supposed to convert nodal data to modal data...
  // ...but I don't need it.  So, just copy the data over.

  int ic = *cellID;
  *nPtsOut = *nPtsIn;

  // Copy uIn to uOut
  uOut = uIn;

  // Get starting offset for global solution-point index
  *iStart = 0;
  for (int i=0; i<ic; i++)
    *iStart += eles[i].getNSpts();
}