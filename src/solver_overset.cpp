/*!
 * \file solver_overset.cpp
 * \brief Overset-related methods for solver class
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

#include "solver.hpp"

/* ---- My New Overset Grid Functions ---- */

void solver::oversetFieldInterp(void)
{
#ifndef _NO_MPI
  if (params->motion) return;  // For moving problems, projection done in moveMesh()

  // Use field interpolation rather than boundary interpolation
  OComm->performProjection_static(eles,Geo->eleMap);
#endif
}

void solver::oversetInterp(void)
{
#ifndef _NO_MPI
  if (params->oversetMethod == 2) return;

  OComm->exchangeOversetData(eles,opers,Geo->eleMap);
#endif
}

void solver::oversetInterp_gradient(void)
{
#ifndef _NO_MPI
  if (params->oversetMethod == 2) return;

  OComm->exchangeOversetGradient(eles,opers,Geo->eleMap);
#endif
}

void solver::setupOverset(void)
{
#ifndef _NO_MPI
  if (gridRank == 0) cout << "Solver: Grid " << gridID << ": Setting up overset connectivity" << endl;

  if (Geo->nDims == 3) {
    OComm = make_shared<overComm>();

    OComm->setup(params,nGrids,gridID,gridRank,nprocPerGrid,Geo->gridIdList);

    OComm->tg = Geo->tg;

    if (params->oversetMethod != 2)
      OComm->matchOversetPoints(eles,overFaces,Geo->eleMap);
  }
  else {
    OComm = Geo->OComm;

    if (params->oversetMethod != 2)
      OComm->matchOversetPoints(eles,overFaces,Geo->eleMap,Geo->minPt,Geo->maxPt);
  }
#endif
}

vector<double> solver::integrateErrorOverset(void)
{
#ifndef _NO_MPI
  return OComm->integrateErrOverset(eles,opers,Geo->iblankCell,Geo->eleMap,params->quadOrder);
#endif
}

/* ---- Basic Tioga-Based Overset-Grid Functions ---- */

void solver::setupOversetData(void)
{
  // Allocate storage for global solution vector (for use with Tioga)
  // and initialize to 0
  int nSptsTotal = 0;
  for (uint i=0; i<eles.size(); i++) nSptsTotal += eles[i]->getNSpts();
  U_spts.assign(nSptsTotal*params->nFields,0);
}

void solver::setGlobalSolutionArray(void)
{
  int ind = 0;
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->getUSpts(&U_spts[ind]);
    ind += eles[i]->getNSpts() * params->nFields;
  }
}

void solver::updateElesSolutionArrays(void)
{
  int ind = 0;
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->setUSpts(&U_spts[ind]);
    ind += eles[i]->getNSpts() * params->nFields;
  }
}

void solver::callDataUpdateTIOGA(void)
{
#ifndef _NO_MPI
  tg->dataUpdate_highorder(params->nFields,U_spts.data(),0);
#endif
}

/* ---- Callback Functions for TIOGA Overset Grid Library ---- */

void solver::getNodesPerCell(int* cellID, int* nNodes)
{
  (*nNodes) = eles[*cellID]->getNSpts();
}

void solver::getReceptorNodes(int* cellID, int* nNodes, double* posNodes)
{
  if (*cellID >= eles.size()) cout << "Invalid cellID!  cellID = " << *cellID << endl;
  eles[*cellID]->getPosSpts(posNodes);
}

void solver::donorInclusionTest(int* cellID, double* xyz, int* passFlag, double* rst)
{
  // Determine if point is in cell: [x,y,z] should all lie between [-1,1]
  point refPt;
  (*passFlag) = eles[*cellID]->getRefLocNelderMead(point(xyz),refPt);

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
    iStart += eles[i]->getNSpts();

  // Put the global indices of ele's spts into inode array
  (*nWeights) = eles[ic]->getNSpts();
  for (int i=0; i<(*nWeights); i++)
    iNode[i] = iStart + i;

  opers[eles[ic]->eType][eles[ic]->order].getBasisValues(rst,weights);
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
    *iStart += eles[i]->getNSpts();
}
