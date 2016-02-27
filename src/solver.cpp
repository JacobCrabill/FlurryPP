/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include <sstream>
#include <omp.h>

class intFace;
class boundFace;

#include "cblas.h"

#include "input.hpp"
#include "geo.hpp"
#include "intFace.hpp"
#include "boundFace.hpp"

solver::solver()
{

}

solver::~solver()
{

}

void solver::setup(input *params, int order, geo *_Geo)
{
  this->params = params;

  if (_Geo == NULL) {
    Geo = new geo;
    Geo->setup(params);
  }
  else {
    Geo = _Geo;
  }

#ifndef _NO_MPI
  this->tg = Geo->tg; // Geo will have initialized this already if needed
#endif

  params->time = 0.;
  this->order = order;

  nDims = params->nDims;
  nFields = params->nFields;
  nMpts = (nDims==2) ? 4 : 8;
  nRKSteps = params->nRKSteps;

  /* Setup the FR elements & faces which will be computed on */
  Geo->setupElesFaces(params,eles,faces,mpiFaces,overFaces);

  nEles = eles.size();

  nGrids = Geo->nGrids;
  gridID = Geo->gridID;
  gridRank = Geo->gridRank;
  nprocPerGrid = Geo->nProcGrid;

  /* Setup the FR operators for computation */
  setupOperators();

  nSpts = opers[order].nSpts;
  nFpts = opers[order].nFpts;

  setupArrays();

  setupGeometry();

  setupElesFaces();

  if (params->meshType == OVERSET_MESH)
    setupOverset();

#ifndef _NO_MPI
  finishMpiSetup();
#endif
}

void solver::setupArrays(void)
{
  U_spts.setup(nSpts, nEles, nFields);
  U_fpts.setup(nFpts, nEles, nFields);

  F_spts.setup(nDims, nSpts, nEles, nFields);
  F_fpts.setup(nDims, nFpts, nEles, nFields);

  if (params->viscous || params->motion)
  {
    dU_spts.setup(nDims, nSpts, nEles, nFields);
    dU_fpts.setup(nDims, nFpts, nEles, nFields);
  }

  if (params->viscous)
  {
    dUc_fpts.setup(nFpts, nEles, nFields);
  }

  dF_spts.setup(nDims, nDims);
  for (auto &mat:dF_spts.data)
    mat.setup(nSpts, nEles, nFields);

  U0.setup(nSpts, nEles, nFields);
  U_mpts.setup(nMpts, nEles, nFields);

  disFn_fpts.setup(nFpts, nEles, nFields);
  Fn_fpts.setup(nFpts, nEles, nFields);

  divF_spts.resize(nRKSteps);
  for (auto &divF:divF_spts)
    divF.setup(nSpts, nEles, nFields);

  /* Multigrid Variables */
  if (params->PMG)
  {
    sol_spts.setup(nSpts, nEles, nFields);
    corr_spts.setup(nSpts, nEles, nFields);
    src_spts.setup(nSpts, nEles, nFields);
  }

  tempVars_spts.setup(nSpts, nEles, nFields);
  tempVars_fpts.setup(nFpts, nEles, nFields);
}

void solver::setupGeometry(void)
{
  nNodes = Geo->nNodesPerCell;

  if (nDims == 2) {
    nPpts = (order+3)*(order+3);
  } else {
    nPpts = (order+3)*(order+3)*(order+3);
  }

  shape_spts.setup(nSpts,nNodes);
  shape_fpts.setup(nFpts,nNodes);

  dshape_spts.setup(nSpts,nNodes,nDims);
  dshape_fpts.setup(nFpts,nNodes,nDims);

  pos_spts.setup(nSpts, nEles, nDims);
  pos_fpts.setup(nFpts, nEles, nDims);
  pos_ppts.setup(nPpts, nEles, nDims);

  tNorm_fpts.setup(nFpts,nDims);

  Jac_spts.setup(nSpts, nEles, nDims, nDims);
  Jac_fpts.setup(nFpts, nEles, nDims, nDims);
  JGinv_spts.setup(nSpts, nEles, nDims, nDims);
  JGinv_fpts.setup(nFpts, nEles, nDims, nDims);
  detJac_spts.setup(nSpts,nEles);
  detJac_fpts.setup(nFpts,nEles);
  dA_fpts.setup(nFpts, nEles);
  norm_fpts.setup(nFpts, nEles, nDims);


  if (params->motion)
  {
    gridV_spts.setup(nSpts, nEles, nDims);
    gridV_fpts.setup(nFpts, nEles, nDims);
    gridV_mpts.setup(nNodes, nEles, nDims);
    gridV_ppts.setup(nPpts, nEles, nDims);
  }

  if (nDims == 2)
  {
    loc_spts = getLocSpts(QUAD,order,params->sptsTypeQuad);
    loc_fpts = getLocFpts(QUAD,order,params->sptsTypeQuad);

    for (uint spt = 0; spt < nSpts; spt++)
    {
      shape_quad(loc_spts[spt], &shape_spts(spt,0),nNodes);
      dshape_quad(loc_spts[spt], &dshape_spts(spt,0,0),nNodes);
    }

    // Setting unit normal vector in the parent domain
    for (uint fpt = 0; fpt < nFpts; fpt++)
    {
      shape_quad(loc_fpts[fpt], &shape_fpts(fpt,0),nNodes);
      dshape_quad(loc_fpts[fpt], &dshape_fpts(fpt,0,0),nNodes);

      uint iFace = floor(fpt / (order+1));
      // Face ordering for quads: Bottom, Right, Top, Left
      switch(iFace) {
        case 0:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = -1;
          break;
        case 1:
          tNorm_fpts(fpt,0) = 1;
          tNorm_fpts(fpt,1) = 0;
          break;
        case 2:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = 1;
          break;
        case 3:
          tNorm_fpts(fpt,0) = -1;
          tNorm_fpts(fpt,1) = 0;
          break;
      }
    }
  }
  else
  {
    loc_spts = getLocSpts(HEX,order,params->sptsTypeQuad);
    loc_fpts = getLocFpts(HEX,order,params->sptsTypeQuad);

    for (uint spt = 0; spt < nSpts; spt++)
    {
      shape_hex(loc_spts[spt], &shape_spts(spt,0), nNodes);
      dshape_hex(loc_spts[spt], &dshape_spts(spt,0,0), nNodes);
    }

    for (uint fpt = 0; fpt < nFpts; fpt++)
    {
      shape_hex(loc_fpts[fpt], &shape_fpts(fpt,0), nNodes);
      dshape_hex(loc_fpts[fpt], &dshape_fpts(fpt,0,0), nNodes);

      // Setting unit normal vector in the parent domain
      uint iFace = floor(fpt / ((order+1)*(order+1)));
      switch(iFace) {
        case 0:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = 0;
          tNorm_fpts(fpt,2) = -1;
          break;
        case 1:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = 0;
          tNorm_fpts(fpt,2) = 1;
          break;
        case 2:
          tNorm_fpts(fpt,0) = -1;
          tNorm_fpts(fpt,1) = 0;
          tNorm_fpts(fpt,2) = 0;
          break;
        case 3:
          tNorm_fpts(fpt,0) = 1;
          tNorm_fpts(fpt,1) = 0;
          tNorm_fpts(fpt,2) = 0;
          break;
        case 4:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = -1;
          tNorm_fpts(fpt,2) = 0;
          break;
        case 5:
          tNorm_fpts(fpt,0) = 0;
          tNorm_fpts(fpt,1) = 1;
          tNorm_fpts(fpt,2) = 0;
          break;
      }
    }
  }

  setPosSptsFpts();
}

void solver::update(bool PMG_Source)
{
  params->iter++;

  /* Intermediate residuals for Runge-Kutta time integration */

  if (params->dtType != 0) calcDt();

  for (int step=0; step<nRKSteps-1; step++) {
    params->rkTime = params->time + params->RKa[step]*params->dt;

    moveMesh(step);

    if (step == 0) copyUspts_U0(); // Store starting values for RK method

    calcResidual(step);

    timeStepA(step, PMG_Source);
  }

  /* Final Runge-Kutta time advancement step */

  params->rkTime = params->time + params->RKa[nRKSteps-1]*params->dt;

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  // Reset solution to initial-stage values
  if (nRKSteps>1)
    copyU0_Uspts();

  for (int step=0; step<nRKSteps; step++)
    timeStepB(step, PMG_Source);

  params->time += params->dt;
}


void solver::calcResidual(int step)
{
  if (params->meshType == OVERSET_MESH && params->oversetMethod == 2) {
    oversetFieldInterp();
  }

  if(params->scFlag == 1) {
    shockCapture();
  }

  extrapolateU();

  /* --- Polynomial-Squeezing stabilization procedure --- */
  if (params->squeeze) {

    calcAvgSolution();

    checkEntropy();

  }

  if (params->viscous || params->motion) {

    calcGradU_spts();

  }

#ifndef _NO_MPI
  doCommunication();
#endif

  calcInviscidFlux_spts();

  calcInviscidFlux_faces();

#ifndef _NO_MPI
  calcInviscidFlux_mpi();
#endif

  if (params->meshType == OVERSET_MESH) {

    oversetInterp();

    calcInviscidFlux_overset();

  }

  if (params->viscous) {

    correctGradU();

    extrapolateGradU();

#ifndef _NO_MPI
    doCommunicationGrad();
#endif

    calcViscousFlux_spts();

    calcViscousFlux_faces();

#ifndef _NO_MPI
    calcViscousFlux_mpi();
#endif

    if (params->meshType == OVERSET_MESH) {
      oversetInterp_gradient();

      calcViscousFlux_overset();
    }
  }

  extrapolateNormalFlux();

  calcFluxDivergence(step);

  correctDivFlux(step);
}

void solver::calcDt(void)
{
  double dt = INFINITY;

#pragma omp parallel for reduction(min:dt)
  for (uint i=0; i<eles.size(); i++) {
    dt = min(dt, eles[i]->calcDt());
  }

#ifndef _NO_MPI
  double dtTmp = dt;
  MPI_Allreduce(&dtTmp, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  params->dt = dt;
}

void solver::timeStepA(int step, bool PMG_Source)
{
  //! TODO
//  if (params->meshType==OVERSET_MESH && params->oversetMethod==2) {
//    for (uint i=0; i<eles.size(); i++) {
//      if (Geo->iblankCell[eles[i]->ID] == NORMAL)

  if (PMG_Source)
  {
    /* --- Include PMG Source Term --- */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        for (uint k = 0; k < nFields; k++) {
          if (params->dtType != 2)
            U_spts(spt,e,k) = U0(spt,e,k) - params->RKa[step+1] * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt,e,k) = U0(spt,e,k) - params->RKa[step+1] * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * eles[e]->dt;
        }
      }
    }
  }
  else
  {
    /* --- Normal Time Advancement --- */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        for (uint k = 0; k < nFields; k++) {
          if (params->dtType != 2)
            U_spts(spt, e, k) = U0(spt,e,k) - params->RKa[step+1] * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt, e, k) = U0(spt,e,k) - params->RKa[step+1] * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * eles[e]->dt;
        }
      }
    }
  }
}

void solver::timeStepB(int step, bool PMG_Source)
{
  if (PMG_Source)
  {
    /* --- Include PMG Source Term --- */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        for (uint k = 0; k < nFields; k++) {
          if (params->dtType != 2)
            U_spts(spt,e,k) -= params->RKb[step] * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt,e,k) -= params->RKb[step] * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * eles[e]->dt;
        }
      }
    }
  }
  else
  {
    /* --- Normal Time Advancement --- */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        for (uint k = 0; k < nFields; k++) {
          if (params->dtType != 2)
            U_spts(spt,e,k) -= params->RKb[step] * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt,e,k) -= params->RKb[step] * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * eles[e]->dt;
        }
      }
    }
  }
}

void solver::copyUspts_U0(void)
{
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++)
      for (uint e = 0; e < nEles; e++)
        for (uint k = 0; k < nFields; k++)
          U0(spt, e, k) = U_spts(spt, e, k);
}

void solver::copyU0_Uspts(void)
{
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++)
      for (uint e = 0; e < nEles; e++)
        for (uint k = 0; k < nFields; k++)
          U_spts(spt, e, k) = U0(spt, e, k);
}

void solver::extrapolateU(void)
{
  int m = nFpts;
  int n = nEles * nFields;
  int k = nSpts;

  auto &A = opers[order].opp_spts_to_fpts(0,0);
  auto &B = U_spts(0,0,0);
  auto &C = U_fpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                    1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif
}

void solver::calcAvgSolution()
{
//#pragma omp parallel for
//  for (uint i=0; i<eles.size(); i++) {
//    opers[order].calcAvgU(eles[i]->U_spts,eles[i]->detJac_spts,eles[i]->Uavg);
//  }
}

bool solver::checkDensity()
{
  bool squeezed = false;
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    bool check = eles[i]->checkDensity();
    squeezed = check|| squeezed;
  }

  return squeezed;
}

void solver::checkEntropy()
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->checkEntropy();
  }
}

void solver::checkEntropyPlot()
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->checkEntropyPlot();
  }
}

void solver::extrapolateUMpts(void)
{
  int m = nMpts;
  int n = nEles * nFields;
  int k = nSpts;

  auto &A = opers[order].opp_spts_to_mpts(0,0);
  auto &B = U_spts(0,0,0);
  auto &C = U_mpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif
}

void solver::extrapolateGridVelMpts(void)
{
//#pragma omp parallel for
//  for (uint i=0; i<eles.size(); i++) {
//    opers[order].applySptsMpts(eles[i]->gridVel_spts,eles[i]->gridVel_mpts);
//  }
}

void solver::extrapolateSMpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[order].applySptsMpts(eles[i]->S_spts,eles[i]->S_mpts);
  }
}

void solver::extrapolateSFpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[order].applySptsFpts(eles[i]->S_spts,eles[i]->S_fpts);
  }
}

void solver::calcInviscidFlux_spts(void)
{
  /*tempF.setup(nFields,nDims);

  double gam1 = params->gamma - 1.;

  if (nDims == 2)
  {
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        double rho = U_spts(spt,e,0);
        double u = U_spts(spt,e,1) / rho;
        double v = U_spts(spt,e,2) / rho;
        double p = gam1*(U_spts(spt,e,3) - 0.5*rho*(u*u + v*v));
        tempF(0,0) =  U_spts(spt,e,1);       tempF(0,1) =  U_spts(spt,e,2);
        tempF(1,0) =  U_spts(spt,e,1)*u+p;   tempF(1,1) =  U_spts(spt,e,1)*v;
        tempF(2,0) =  U_spts(spt,e,2)*u;     tempF(2,1) =  U_spts(spt,e,2)*v+p;
        tempF(3,0) = (U_spts(spt,e,3)+p)*u;  tempF(3,1) = (U_spts(spt,e,3)+p)*v;

        / * --- Transform back to reference domain --- * /
        for (uint dim1 = 0; dim1 < nDims; dim1++) {
          for (uint k = 0; k < nFields; k++) {
            F_spts(dim1,spt,e,k) = 0.;
            for (uint dim2 = 0; dim2 < nDims; dim2++) {
              F_spts(dim1, spt,e, k) += JGinv_spts(spt,e,dim1,dim2)*tempF(k,dim2);
            }
          }
        }
      }
    }
  }*/
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->calcInviscidFlux_spts();
  }
}

void solver::doCommunication()
{
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->communicate();
  }
}

void solver::doCommunicationGrad()
{
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->communicateGrad();
  }
}

void solver::calcInviscidFlux_faces()
{
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->calcInviscidFlux();
  }
}

void solver::calcInviscidFlux_mpi()
{
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->calcInviscidFlux();
  }
}

void solver::calcInviscidFlux_overset()
{
  if (params->oversetMethod == 2) return;

#pragma omp parallel for
  for (uint i=0; i<overFaces.size(); i++) {
    overFaces[i]->calcInviscidFlux();
  }
}

void solver::calcViscousFlux_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->calcViscousFlux_spts();
  }
}

void solver::calcViscousFlux_faces()
{
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->calcViscousFlux();
  }
}

void solver::calcViscousFlux_mpi()
{
  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->calcViscousFlux();
  }
}

void solver::calcViscousFlux_overset()
{
  if (params->oversetMethod == 2) return;

#pragma omp parallel for
  for (uint i=0; i<overFaces.size(); i++) {
    overFaces[i]->calcViscousFlux();
  }
}

void solver::calcGradF_spts(void)
{
  int m = nSpts;
  int n = nEles * nFields;
  int k = nSpts;

  for (uint dim1=0; dim1<nDims; dim1++) {
    for (uint dim2=0; dim2<nDims; dim2++) {
      auto &A = opers[order].opp_grad_spts[dim2](0,0);
      auto &B = F_spts(dim1,0,0,0);
      auto &C = dF_spts(dim2,dim1)(0,0,0);
#ifdef _OMP
      omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 0.0, &C, n);
#else
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 0.0, &C, n);
#endif
    }
  }
}

void solver::transformGradF_spts(int step)
{
  divF_spts[step].initializeToValue(0);

  if (nDims == 2)
  {
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        double A = gridV_spts(spt,e,1)*Jac_spts(spt,e,0,1) - gridV_spts(spt,e,0)*Jac_spts(spt,e,1,1);
        double B = gridV_spts(spt,e,0)*Jac_spts(spt,e,1,0) - gridV_spts(spt,e,1)*Jac_spts(spt,e,0,0);
        for (uint k = 0; k < nFields; k++) {
          dF_spts(0,0)(spt,e,k) =  dF_spts(0,0)(spt,e,k)*Jac_spts(spt,e,1,1) - dF_spts(0,1)(spt,e,k)*Jac_spts(spt,e,0,1) + dU_spts(0,spt,e,k)*A;
          dF_spts(1,1)(spt,e,k) = -dF_spts(1,0)(spt,e,k)*Jac_spts(spt,e,1,0) + dF_spts(1,1)(spt,e,k)*Jac_spts(spt,e,0,0) + dU_spts(1,spt,e,k)*B;
          divF_spts[step](spt,e,k) = dF_spts(0,0)(spt,e,k) + dF_spts(1,1)(spt,e,k);
        }
      }
    }
  }
  else
  {
    matrix<double> Jacobian(4,4);
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        Jacobian(3,3) = 1;
        for (uint i = 0; i < 3; i++) {
          for (uint j = 0; j < 3; j++)
            Jacobian(i,j) = Jac_spts(spt,e,i,j);
          Jacobian(i,3) = gridV_spts(spt,e,i);
        }
        matrix<double> S = Jacobian.adjoint();

        for (uint dim1 = 0; dim1 < 3; dim1++)
          for (uint dim2 = 0; dim2 < 3; dim2++)
            for (uint k = 0; k < nFields; k++)
              divF_spts[step](spt,e,k) += dF_spts(dim2,dim1)(spt,e,k)*S(dim2,dim1);

        for (uint dim = 0; dim < 3; dim++)
          for (uint k = 0; k < nFields; k++)
            divF_spts[step](spt,e,k) += dU_spts(dim,spt,e,k)*S(dim,3);
      }
    }
  }
}

void solver::calcFluxDivergence(int step)
{
  if (params->motion) {

    /* Use non-conservation-form chain-rule transformation
     * (See AIAA paper 2013-0998 by Liang, Miyaji and Zhang) */

    calcGradF_spts();

    transformGradF_spts(step);

  } else {

    /* Standard conservative form */
    calcDivF_spts(step);

  }
}

void solver::calcDivF_spts(int step)
{
  int m = nSpts;
  int n = nEles * nFields;
  int k = nSpts;

  auto &C = divF_spts[step](0,0,0);

  auto &A = opers[order].opp_grad_spts[0](0,0);
  auto &B = F_spts(0,0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif

  for (uint dim = 1; dim < nDims; dim++) {
    auto &A = opers[order].opp_grad_spts[dim](0,0);
    auto &B = F_spts(dim,0,0,0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 1.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 1.0, &C, n);
#endif
  }
}

void solver::extrapolateNormalFlux(void)
{
  if (params->motion)
  {
    /* Extrapolate physical normal flux */

    int m = nFpts;
    int n = nEles * nFields;
    int k = nSpts;

    auto &A = opers[order].opp_spts_to_fpts(0,0);

    auto &B = F_spts(0, 0, 0, 0);
    auto &C = disFn_fpts(0, 0, 0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#endif

#pragma omp parallel for collapse(3)
    for (uint fpt = 0; fpt < nFpts; fpt++)
      for (uint e = 0; e < nEles; e++)
        for (uint k = 0; k < nFields; k++)
          disFn_fpts(fpt,e,k) *= eles[e]->norm_fpts(fpt,0) * dA_fpts(fpt,e);

    for (uint dim = 1; dim < nDims; dim++) {
      auto &B = F_spts(dim, 0, 0, 0);
      auto &C = tempVars_fpts(0, 0, 0);
#ifdef _OMP
      omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 0.0, &C, n);
#else
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 0.0, &C, n);
#endif

#pragma omp parallel for collapse(3)
      for (uint fpt = 0; fpt < nFpts; fpt++)
        for (uint e = 0; e < nEles; e++)
          for (uint k = 0; k < nFields; k++)
            disFn_fpts(fpt,e,k) += tempVars_fpts(fpt,e,k) * eles[e]->norm_fpts(fpt,dim) * dA_fpts(fpt,e);
    }
  }
  else
  {
    /* Extrapolate transformed normal flux */

    int m = nFpts;
    int n = nEles * nFields;
    int k = nSpts;

    auto &C = disFn_fpts(0, 0, 0);

    auto &A = opers[order].opp_extrapolateFn[0](0,0);
    auto &B = F_spts(0, 0, 0, 0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#endif

    for (uint dim = 1; dim < nDims; dim++) {
      auto &A = opers[order].opp_extrapolateFn[dim](0,0);
      auto &B = F_spts(dim, 0, 0, 0);
#ifdef _OMP
      omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 1.0, &C, n);
#else
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                  1.0, &A, k, &B, n, 1.0, &C, n);
#endif
    }
  }
}

void solver::correctDivFlux(int step)
{
#pragma omp parallel for collapse(3)
  for (uint fpt = 0; fpt < nFpts; fpt++)
    for (uint e = 0; e < nEles; e++)
      for (uint k = 0; k < nFields; k++)
        disFn_fpts(fpt, e, k) = Fn_fpts(fpt, e, k) - disFn_fpts(fpt, e, k);

  int m = nSpts;
  int n = nEles * nFields;
  int k = nFpts;

  auto &A = opers[order].opp_correction(0,0);
  auto &B = disFn_fpts(0,0,0);
  auto &C = divF_spts[step](0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 1.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 1.0, &C, n);
#endif
}

void solver::calcGradU_spts(void)
{
  int m = nSpts;
  int n = nEles * nFields;
  int k = nSpts;

  for (uint dim1=0; dim1<nDims; dim1++) {
    auto &A = opers[order].opp_grad_spts[dim1](0,0);
    auto &B = U_spts(0,0,0);
    auto &C = dU_spts(dim1,0,0,0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#endif
  }
}

void solver::correctGradU(void)
{
  /* Apply correction to solution gradient in reference space */

  int m = nSpts;
  int n = nFpts;
  int k = nEles * nFields;

  auto &B = dUc_fpts(0,0,0);

  auto &A = opers[order].opp_correctU[0](0,0);
  auto &C = dU_spts(0, 0, 0, 0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif

  for (uint dim = 0; dim < nDims; dim++) {
    auto &A = opers[order].opp_correctU[dim](0,0);
    auto &C = dU_spts(dim, 0, 0, 0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 1.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 1.0, &C, n);
#endif
  }

  /* Transform back to physical space */

  if (nDims == 2)
  {
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        double invDet = 1./detJac_spts(spt,e);
        for (uint k = 0; k < nFields; k++) {
          double ur = dU_spts(0,spt,e,k);
          double us = dU_spts(1,spt,e,k);
          dU_spts(0,spt,e,k) = invDet * (ur*JGinv_spts(spt,e,0,0) + us*JGinv_spts(spt,e,1,0));
          dU_spts(1,spt,e,k) = invDet * (ur*JGinv_spts(spt,e,0,1) + us*JGinv_spts(spt,e,1,1));
        }
      }
    }
  }
  else
  {
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        double invDet = 1./detJac_spts(spt,e);
        for (uint k = 0; k < nFields; k++) {
          double ur = dU_spts(0,spt,e,k);
          double us = dU_spts(1,spt,e,k);
          double ut = dU_spts(2,spt,e,k);
          dU_spts(0,spt,e,k) = invDet * (ur*JGinv_spts(spt,e,0,0) + us*JGinv_spts(spt,e,1,0) + ut*JGinv_spts(spt,e,2,0));
          dU_spts(1,spt,e,k) = invDet * (ur*JGinv_spts(spt,e,0,1) + us*JGinv_spts(spt,e,1,1) + ut*JGinv_spts(spt,e,2,1));
          dU_spts(2,spt,e,k) = invDet * (ur*JGinv_spts(spt,e,0,2) + us*JGinv_spts(spt,e,1,2) + ut*JGinv_spts(spt,e,2,2));
        }
      }
    }
  }
}

void solver::extrapolateGradU()
{
  int m = nFpts;
  int n = nEles * nFields;
  int k = nSpts;

  auto &A = opers[order].opp_spts_to_fpts(0,0);

  for (uint dim=0; dim<nDims; dim++) {
    auto &B = dU_spts(dim,0,0,0);
    auto &C = dU_fpts(dim,0,0,0);
#ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
#endif
  }
}

void solver::calcEntropyErr_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->calcEntropyErr_spts();
  }
}

void solver::moveMesh(int step)
{
  if (!params->motion) return;

  if (params->meshType == OVERSET_MESH) {
    if (step == 0) {
      Geo->setIterIblanks();
      if (params->RKa[nRKSteps-1]!=1.)
        Geo->updateADT();
      Geo->processBlanks(eles,faces,mpiFaces,overFaces);
      Geo->processUnblanks(eles,faces,mpiFaces,overFaces,this);
      OComm->matchUnblankCells(eles,Geo->unblankCells,Geo->eleMap,params->quadOrder);
      OComm->performGalerkinProjection(eles,opers,Geo->eleMap,order);
    }

    if (params->oversetMethod == 2) {
      if (params->projection) {
        OComm->matchUnblankCells(eles,Geo->fringeCells,Geo->eleMap,params->quadOrder);
        OComm->performGalerkinProjection(eles,opers,Geo->eleMap,order);
      } else {
        if (params->nDims==2)
          getBoundingBox(Geo->xv,Geo->minPt,Geo->maxPt);
        OComm->setupFringeCellPoints(eles,Geo->fringeCells,Geo->eleMap);
        OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
        OComm->exchangeOversetData(eles,opers,Geo->eleMap);
        OComm->transferEleData(eles,Geo->fringeCells,Geo->eleMap);
      }
    }

    if ( !(step==0 && params->RKa[step]==0) )
      Geo->moveMesh(params->RKa[step]);

    updatePosSptsFpts();

    updateGridVSptsFpts();

    if (params->motion != 4) {
#pragma omp parallel for
      for (uint i=0; i<eles.size(); i++)
        eles[i]->calcTransforms(true);
    }

    Geo->updateADT();

    if (params->oversetMethod != 2) {
      if (params->nDims==2)
        getBoundingBox(Geo->xv,Geo->minPt,Geo->maxPt);
      OComm->setupOverFacePoints(overFaces);
      OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
    }
  } else {
    Geo->moveMesh(params->RKa[step]);

    updatePosSptsFpts();

    updateGridVSptsFpts();

    updateTransforms();
  }

}

void solver::setPosSptsFpts(void)
{
  nodes.setup(nNodes, nEles, nDims);

  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        nodes(npt, e, dim) = Geo->xv[Geo->c2v(e,npt)][dim];

  int ms = nSpts;
  int mf = nFpts;
  int k = nNodes;
  int n = nEles * nDims;

  auto &B = nodes(0,0,0);

  auto &As = shape_spts(0,0);
  auto &Cs = pos_spts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#endif

  auto &Af = shape_fpts(0,0);
  auto &Cf = pos_fpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#endif

  /* Initialize storage of moving node positions */
  if (params->motion)
  {
    nodesRK = nodes;
  }
}

void solver::updatePosSptsFpts(void)
{
#pragma omp parallel for collapse(3)
  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        nodesRK(npt, e, dim) = Geo->xv[Geo->c2v(e,npt)][dim];

  int ms = nSpts;
  int mf = nFpts;
  int k = nNodes;
  int n = nEles * nDims;

  auto &B = nodesRK(0,0,0);

  auto &As = shape_spts(0,0);
  auto &Cs = pos_spts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#endif

  auto &Af = shape_fpts(0,0);
  auto &Cf = pos_fpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#endif
}

void solver::updateGridVSptsFpts(void)
{
#pragma omp parallel for collapse(3)
  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        gridV_mpts(npt, e, dim) = Geo->gridVel(Geo->c2v(e,npt),dim);

  int ms = nSpts;
  int mf = nFpts;
  int k = nNodes;
  int n = nEles * nDims;

  auto &B = gridV_mpts(0,0,0);

  auto &As = shape_spts(0,0);
  auto &Cs = gridV_spts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, n, k,
              1.0, &As, k, &B, n, 0.0, &Cs, n);
#endif

  auto &Af = shape_fpts(0,0);
  auto &Cf = gridV_fpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, n, k,
              1.0, &Af, k, &B, n, 0.0, &Cf, n);
#endif
}

void solver::calcTransforms(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++)
    eles[i]->calcTransforms(false);
}

void solver::updateTransforms(void)
{
  if (params->motion != 4) {
#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++)
      eles[i]->calcTransforms(true);
  }
}

vector<double> solver::computeWallForce(void)
{
  vector<double> force = {0,0,0,0,0,0};

  for (uint i=0; i<faces.size(); i++) {
    auto fTmp = faces[i]->computeWallForce();

    for (int j=0; j<6; j++)
      force[j] += fTmp[j];
  }

  return force;
}

vector<double> solver::computeMassFlux(void)
{
  vector<double> flux(params->nFields);

  for (uint i=0; i<faces.size(); i++) {
    auto fTmp = faces[i]->computeMassFlux();

    for (int j=0; j<params->nFields; j++)
      flux[j] += fTmp[j];
  }

#ifndef _NO_MPI
    auto fTmp = flux;
    MPI_Reduce(fTmp.data(), flux.data(), params->nFields, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

  return flux;
}

void solver::setupOperators()
{
  if (params->rank==0) cout << "Solver: Setting up FR operators" << endl;

  if (params->nDims == 2)
    opers[order].setupOperators(QUAD,order,Geo,params);
  else
    opers[order].setupOperators(HEX,order,Geo,params);
}

void solver::setupElesFaces(void) {

  if (params->rank==0) cout << "Solver: Setting up elements & faces" << endl;

#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i]->setup(params,this,Geo,order);
  }

  // Finish setting up internal & boundary faces
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i]->setupFace();
  }

  // Finish setting up MPI faces

  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->setupFace();
  }

  // Finish setting up overset faces

  for (uint i=0; i<overFaces.size(); i++) {
    overFaces[i]->setupFace();
  }
}

void solver::finishMpiSetup(void)
{
  if (params->rank==0) cout << "Solver: Setting up MPI face communications" << endl;

  for (uint i=0; i<mpiFaces.size(); i++) {
    mpiFaces[i]->finishRightSetup();
  }
}

void solver::readRestartFile(void) {

  ifstream dataFile;
  dataFile.precision(15);

  // Get the file name & open the file
  char fileNameC[256];
  string fileName = params->dataFileName;
#ifndef _NO_MPI
  /* --- All processors read their data from their own .vtu file --- */
  if (params->meshType == OVERSET_MESH)
    sprintf(fileNameC,"%s_%.09d/%s%d_%.09d_%d.vtu",&fileName[0],params->restartIter,&fileName[0],gridID,params->restartIter,gridRank);
  else
    sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],params->restartIter,&fileName[0],params->restartIter,params->rank);
#else
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],params->restartIter);
#endif

  if (params->rank==0) cout << "Solver: Restarting from " << fileNameC << endl;

  dataFile.open(fileNameC);

  if (!dataFile.is_open())
    FatalError("Cannont open restart file.");

  // Read the simulation time from the comment section, then find the start of
  // the UnstructuredData region
  // Also read overset iblank data if applicable
  bool foundTime  = false;
  bool foundIBTag = false;
  bool foundUGTag = false;
  string str;
  stringstream ss;
  vector<double> tmpIblank;
  while (getline(dataFile,str)) {
    ss.str(string("")); ss.clear();
    ss.str(str);
    ss >> str;
    if (str.compare("<!--")==0) {
      ss >> str;
      if (str.compare("TIME")==0) {
        foundTime = true;
        ss >> params->time;
        params->rkTime = params->time;
        if (params->rank == 0)
          cout << "  Restart time = " << params->time << endl;
      } else if (params->meshType == OVERSET_MESH and str.compare("IBLANK_CELL")==0) {
        foundIBTag = true;
        // Read cell Iblank data for overset cases
        Geo->iblankCell.resize(Geo->nEles);
        tmpIblank.assign(Geo->nEles,NORMAL);
        for (int i=0; i<Geo->nEles; i++)
          ss >> tmpIblank[i];
      }
    } else if (str.compare("<UnstructuredGrid>")==0) {
      foundUGTag = true;
      break;
    }
  }

  if (!foundTime)
    cout << "WARNING: Unable to read simulation restart time." << endl;

  if (!foundUGTag)
    FatalError("Cannot find UnstructuredData tag in restart file.");

  if (params->meshType == OVERSET_MESH and !foundIBTag)
    cout << "WARNING: IblankCell data not found in restart file for rank " << params->rank << endl;

  /* -- Set the geometry to the current restart time -- */

  moveMesh(0);

  if (params->meshType == OVERSET_MESH) {
    Geo->unblankCells.clear();
    Geo->blankCells.clear();

    for (int ic=0; ic<Geo->nEles; ic++) {
      Geo->iblankCell[ic] = tmpIblank[ic];
      if (tmpIblank[ic] == HOLE)
        Geo->blankCells.insert(ic);
    }

    Geo->processBlanks(eles,faces,mpiFaces,overFaces);
  }

  // Read restart data & setup all data arrays
  if (params->meshType == OVERSET_MESH && params->oversetMethod == 2) {
    for (int ic=0; ic<eles.size(); ic++) {
      if (tmpIblank[eles[ic]->ID] != NORMAL) continue;
      eles[ic]->restart(dataFile,params,Geo);
    }
  }
  else {
    for (auto& e:eles)
      e->restart(dataFile,params,Geo);
  }

  dataFile.close();

  if (params->rank==0) cout << "Solver: Done reading restart file." << endl;
}

void solver::initializeSolution(bool PMG)
{
  if (params->rank==0) cout << "Solver: Initializing Solution... " << flush;

  if (params->restart && !PMG) {

    readRestartFile();

  } else {

    /* Get the initial grid velocity for wave speed calculations */
    if (params->motion > 0)
      moveMesh(0);

#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++) {
      eles[i]->setInitialCondition();
    }

  }

  if (params->meshType == OVERSET_MESH && params->motion == 0) {
    // Perform initial LGP to setup connectivity / arrays for remainder of computations [Field-interp method]
    if (params->projection) {
      OComm->matchUnblankCells(eles,Geo->fringeCells,Geo->eleMap,params->quadOrder);
      OComm->performGalerkinProjection(eles,opers,Geo->eleMap,order);
    } else {
      if (params->nDims==2)
        getBoundingBox(Geo->xv,Geo->minPt,Geo->maxPt);
      OComm->setupFringeCellPoints(eles,Geo->fringeCells,Geo->eleMap);
      OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
      OComm->exchangeOversetData(eles,opers,Geo->eleMap);
    }
  }

  /* If using CFL-based time-stepping, calc wave speed in each
   * ele for initial dt calculation */
  if (params->dtType != 0) {
    extrapolateU();
#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++) {
      eles[i]->calcWaveSpFpts();
    }
  }

  if (params->rank == 0) cout << "done." << endl;
}

vector<double> solver::integrateError(void)
{
  int quadOrder = params->quadOrder;
  auto wts = getQptWeights(quadOrder,params->nDims);

  /* Interpolate solution to quadrature points */

  auto &qpts = opers[order].loc_qpts;
  int nQpts = qpts.size();

  U_qpts.setup(nQpts, nEles, nFields);
  detJac_qpts.setup(nQpts, nEles);

  int m = nQpts;
  int n = nEles * nFields;
  int k = nSpts;

  auto &A = opers[order].opp_spts_to_qpts(0,0);
  auto &B = U_spts(0,0,0);
  auto &C = U_qpts(0,0,0);

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);

  n = nEles;
  auto &B1 = detJac_spts(0,0);
  auto &C1 = detJac_qpts(0,0);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B1, n, 0.0, &C1, n);

  /* Integrate error over each element */

  vector<double> LpErr(params->nFields);
  for (uint ic = 0; ic < eles.size(); ic++) {
    //if (params->meshType == OVERSET_MESH and Geo->iblankCell[eles[ic]->ID]!=NORMAL) continue;
    for (uint qpt = 0; qpt < nQpts; qpt++) {
      auto tmpErr = calcError(&U_qpts(qpt, ic, 0), eles[ic]->calcPos(qpts[qpt]), params);
      for (uint j = 0; j < params->nFields; j++)
        LpErr[j] += tmpErr[j] * wts[qpt] * detJac_qpts(qpt, ic);
    }
  }

#ifndef _NO_MPI
  vector<double> tmpErr = LpErr;
  MPI_Allreduce(tmpErr.data(), LpErr.data(), params->nFields, MPI_DOUBLE, MPI_SUM, Geo->gridComm);
#endif

  if (params->errorNorm==2)
    for (auto &val:LpErr) val = std::sqrt(std::abs(val));

  return LpErr;
}

// Method for shock capturing
void solver::shockCapture(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
/////    eles[i]->sensor = opers[order].shockCaptureInEle(eles[i]->U_spts,params->threshold);
  }
}
