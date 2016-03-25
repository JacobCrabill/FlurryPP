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
  if (params->projection)
    OComm->performProjection_static(eles,Geo->eleMap,order);
  else {
    OComm->exchangeOversetData(eles,opers,Geo->eleMap);
    OComm->transferEleData(eles,Geo->fringeCells,Geo->eleMap);
  }
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
  }
  else {
    OComm = Geo->OComm;

    getBoundingBox(Geo->xv,Geo->minPt,Geo->maxPt);
  }

  if (params->oversetMethod == 2 && !params->projection) {
    OComm->setupFringeCellPoints(eles,Geo->fringeCells,Geo->eleMap);
  } else if (params->oversetMethod != 2) {
    int nPtsFace = order+1;
    if (nDims==3) nPtsFace *= order+1;
    OComm->setupOverFacePoints(overFaces,nPtsFace);
  }

  OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
#endif
}

vector<double> solver::integrateErrorOverset(void)
{
#ifndef _NO_MPI
  return OComm->integrateErrOverset(eles,opers,Geo->iblankCell,Geo->eleMap,order,params->quadOrder);
#endif
}

/* ---- Basic Tioga-Based Overset-Grid Functions ---- */

//void solver::setupOversetData(void)
//{
//  // Allocate storage for global solution vector (for use with Tioga)
//  // and initialize to 0
//  int nSptsTotal = 0;
//  for (uint i=0; i<eles.size(); i++) nSptsTotal += eles[i]->getNSpts();
//  U_spts.assign(nSptsTotal*params->nFields,0);
//}

//void solver::setGlobalSolutionArray(void)
//{
//  int ind = 0;
//  for (uint i=0; i<eles.size(); i++) {
//    eles[i]->getUSpts(&U_spts[ind]);
//    ind += eles[i]->getNSpts() * params->nFields;
//  }
//}

//void solver::updateElesSolutionArrays(void)
//{
//  int ind = 0;
//  for (uint i=0; i<eles.size(); i++) {
//    eles[i]->setUSpts(&U_spts[ind]);
//    ind += eles[i]->getNSpts() * params->nFields;
//  }
//}

//void solver::callDataUpdateTIOGA(void)
//{
//#ifndef _NO_MPI
//  tg->dataUpdate_highorder(params->nFields,U_spts.data(),0);
//#endif
//}

void solver::insertElement(uint ele_ind)
{
  U_spts.add_dim_1(ele_ind, 0.);
  U_fpts.add_dim_1(ele_ind, 0.);

  F_spts.add_dim_2(ele_ind, 0.);
  F_fpts.add_dim_2(ele_ind, 0.);

  if (params->viscous || params->motion)
  {
    dU_spts.add_dim_2(ele_ind, 0.);
    dU_fpts.add_dim_2(ele_ind, 0.);
  }

  if (params->viscous)
  {
    dUc_fpts.add_dim_1(ele_ind, 0.);
  }

  //dF_spts.setup(nDims, nDims);
  for (auto &mat:dF_spts.data)
    mat.add_dim_1(ele_ind, 0.);

  U0.add_dim_1(ele_ind, 0.);
  U_mpts.add_dim_1(ele_ind, 0.);

  disFn_fpts.add_dim_1(ele_ind, 0.);
  Fn_fpts.add_dim_1(ele_ind, 0.);

  for (auto &divF:divF_spts)
    divF.add_dim_1(ele_ind, 0.);

  /* Multigrid Variables */
  if (params->PMG)
  {
    sol_spts.add_dim_1(ele_ind, 0.);
    corr_spts.add_dim_1(ele_ind, 0.);
    src_spts.add_dim_1(ele_ind, 0.);
  }

  tempVars_spts.add_dim_1(ele_ind, 0.);
  tempVars_fpts.add_dim_1(ele_ind, 0.);

  pos_spts.add_dim_1(ele_ind, 0.);
  pos_fpts.add_dim_1(ele_ind, 0.);
  pos_ppts.add_dim_1(ele_ind, 0.);

  Jac_spts.add_dim_1(ele_ind, 0.);
  Jac_fpts.add_dim_1(ele_ind, 0.);
  JGinv_spts.add_dim_1(ele_ind, 0.);
  JGinv_fpts.add_dim_1(ele_ind, 0.);
  detJac_spts.add_dim_1(ele_ind, 0.);
  detJac_fpts.add_dim_1(ele_ind, 0.);
  dA_fpts.add_dim_1(ele_ind, 0.);
  norm_fpts.add_dim_1(ele_ind, 0.);

  nodes.add_dim_1(ele_ind, 0.);

  if (params->motion)
  {
    gridV_spts.add_dim_1(ele_ind, 0.);
    gridV_fpts.add_dim_1(ele_ind, 0.);
    gridV_mpts.add_dim_1(ele_ind, 0.);
    gridV_ppts.add_dim_1(ele_ind, 0.);

    nodesRK.add_dim_1(ele_ind, 0.);
  }

  nEles++;
}

void solver::removeElement(uint ele_ind)
{
  U_spts.remove_dim_1(ele_ind);
  U_fpts.remove_dim_1(ele_ind);

  F_spts.remove_dim_2(ele_ind);
  F_fpts.remove_dim_2(ele_ind);

  if (params->viscous || params->motion)
  {
    dU_spts.remove_dim_2(ele_ind);
    dU_fpts.remove_dim_2(ele_ind);
  }

  if (params->viscous)
  {
    dUc_fpts.remove_dim_1(ele_ind);
  }

  //dF_spts.setup(nDims, nDims);
  for (auto &mat:dF_spts.data)
    mat.remove_dim_1(ele_ind);

  U0.remove_dim_1(ele_ind);
  U_mpts.remove_dim_1(ele_ind);

  disFn_fpts.remove_dim_1(ele_ind);
  Fn_fpts.remove_dim_1(ele_ind);

  for (auto &divF:divF_spts)
    divF.remove_dim_1(ele_ind);

  /* Multigrid Variables */
  if (params->PMG)
  {
    sol_spts.remove_dim_1(ele_ind);
    corr_spts.remove_dim_1(ele_ind);
    src_spts.remove_dim_1(ele_ind);
  }

  tempVars_spts.remove_dim_1(ele_ind);
  tempVars_fpts.remove_dim_1(ele_ind);

  pos_spts.remove_dim_1(ele_ind);
  pos_fpts.remove_dim_1(ele_ind);
  pos_ppts.remove_dim_1(ele_ind);

  Jac_spts.remove_dim_1(ele_ind);
  Jac_fpts.remove_dim_1(ele_ind);
  JGinv_spts.remove_dim_1(ele_ind);
  JGinv_fpts.remove_dim_1(ele_ind);
  detJac_spts.remove_dim_1(ele_ind);
  detJac_fpts.remove_dim_1(ele_ind);
  dA_fpts.remove_dim_1(ele_ind);
  norm_fpts.remove_dim_1(ele_ind);

  nodes.remove_dim_1(ele_ind);

  if (params->motion)
  {
    gridV_spts.remove_dim_1(ele_ind);
    gridV_fpts.remove_dim_1(ele_ind);
    gridV_mpts.remove_dim_1(ele_ind);
    gridV_ppts.remove_dim_1(ele_ind);

    nodesRK.remove_dim_1(ele_ind);
  }

  nEles--;
}

/* ---- Callback Functions for TIOGA Overset Grid Library ---- */

void solver::getNodesPerCell(int* cellID, int* nNodes)
{
  int ie = Geo->eleMap[*cellID];
  (*nNodes) = eles[ie]->getNSpts();
}

void solver::getReceptorNodes(int* cellID, int* nNodes, double* posNodes)
{
  int ie = Geo->eleMap[*cellID];
  //if (ie >= eles.size()) cout << "Invalid cellID!  cellID = " << *cellID << endl;
  eles[ie]->getPosSpts(posNodes);
}

void solver::donorInclusionTest(int* cellID, double* xyz, int* passFlag, double* rst)
{
  // Determine if point is in cell: [x,y,z] should all lie between [-1,1]
  point refPt = {99., 99., 99.};
  int ie = Geo->eleMap[*cellID];
  if (ie > 0)
    (*passFlag) = eles[ie]->getRefLocNewton(point(xyz),refPt);
  else
    *passFlag = 0;

  rst[0] = refPt.x;
  rst[1] = refPt.y;
  rst[2] = refPt.z;
}

void solver::donorWeights(int* cellID, double* xyz, int* nWeights, int* iNode, double* weights, double* rst, int* fracSize)
{
  int ic = Geo->eleMap[*cellID];

  // Get starting offset for global solution-point index
  int iStart = 0;
  for (int i=0; i<ic; i++)
    iStart += eles[i]->getNSpts();

  // Put the global indices of ele's spts into inode array
  (*nWeights) = eles[ic]->getNSpts();
  for (int i=0; i<(*nWeights); i++)
    iNode[i] = iStart + i;

  opers[eles[ic]->order].getBasisValues(rst,weights);
}

void solver::convertToModal(int* cellID, int* nPtsIn, double* uIn, int* nPtsOut, int* iStart, double* uOut)
{
  // This is apparently supposed to convert nodal data to modal data...
  // ...but I don't need it.  So, just copy the data over.

  int ic = Geo->eleMap[*cellID];
  *nPtsOut = *nPtsIn;

  // Copy uIn to uOut
  uOut = uIn;

  // Get starting offset for global solution-point index
  *iStart = 0;
  for (int i=0; i<ic; i++)
    *iStart += eles[i]->getNSpts();
}
