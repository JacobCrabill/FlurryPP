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

#ifdef _MKL_BLAS
#include "mkl_cblas.h"
#else
#include "cblas.h"
#endif

#include "flux.hpp"
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

void solver::setup(input *params, int _order, geo *_Geo)
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
  order = _order;

  nDims = params->nDims;
  nFields = params->nFields;
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
  nPpts = opers[order].nPpts;
  nMpts = nPpts - nSpts - nFpts;
  nCpts = nSpts;

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
  U_mpts.setup(nMpts, nEles, nFields);
  V_ppts.setup(nPpts, nEles, nFields);
  V_spts.setup(nSpts,nEles,nFields); /// TEMP / DEBUGGING

  F_spts.setup(nDims, nSpts, nEles, nFields);
  F_fpts.setup(nDims, nFpts, nEles, nFields);

  if (params->viscous || params->motion)
  {
    dU_spts.setup(nDims, nSpts, nEles, nFields);
    dU_fpts.setup(nDims, nFpts, nEles, nFields);
    tempDU.setup(nDims,nFields);
  }

  if (params->viscous)
  {
    dUc_fpts.setup(nFpts, nEles, nFields);
  }

  dF_spts.setup(nDims, nDims);
  for (auto &mat:dF_spts.data)
    mat.setup(nSpts, nEles, nFields);

  U0.setup(nSpts, nEles, nFields);

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
  nNodes = getMax(Geo->c2nv);

  shape_spts.setup(nSpts,nNodes);
  shape_fpts.setup(nFpts,nNodes);
  shape_ppts.setup(nPpts,nNodes);
  shape_cpts.setup(nCpts,nNodes);

  dshape_spts.setup(nDims,nSpts,nNodes);
  dshape_fpts.setup(nDims,nFpts,nNodes);
  dshape_cpts.setup(nDims,nCpts,nNodes);

  pos_spts.setup(nSpts, nEles, nDims);
  pos_fpts.setup(nFpts, nEles, nDims);
  pos_ppts.setup(nPpts, nEles, nDims);
  pos_cpts.setup(nCpts, nEles, nDims);

  tNorm_fpts.setup(nFpts,nDims);

  Jac_spts.setup(nSpts, nEles, nDims, nDims);
  Jac_fpts.setup(nFpts, nEles, nDims, nDims);
  JGinv_spts.setup(nSpts, nEles, nDims, nDims);
  JGinv_fpts.setup(nFpts, nEles, nDims, nDims);
  detJac_spts.setup(nSpts,nEles);
  detJac_fpts.setup(nFpts,nEles);
  dA_fpts.setup(nFpts, nEles);
  norm_fpts.setup(nFpts, nEles, nDims);

  nodes.setup(nNodes, nEles, nDims);

  if (params->motion)
  {
    gridV_spts.setup(nSpts, nEles, nDims);
    gridV_fpts.setup(nFpts, nEles, nDims);
    gridV_mpts.setup(nNodes, nEles, nDims);
    gridV_ppts.setup(nPpts, nEles, nDims);

    nodesRK.setup(nNodes, nEles, nDims);
  }

  if (nDims == 2)
  {
    loc_spts = getLocSpts(QUAD,order,params->sptsTypeQuad);
    loc_fpts = getLocFpts(QUAD,order,params->sptsTypeQuad);
    loc_ppts = getLocPpts(QUAD,order,params->sptsTypeQuad);
    loc_cpts = getLocSpts(QUAD,order,string("Equidistant"));

    dshape_quad(loc_spts, dshape_spts, nNodes);
    dshape_quad(loc_fpts, dshape_fpts, nNodes);
    dshape_quad(loc_cpts, dshape_cpts, nNodes);

    for (uint ppt = 0; ppt < nPpts; ppt++)
      shape_quad(loc_ppts[ppt], &shape_ppts(ppt,0), nNodes);

    for (uint cpt = 0; cpt < nCpts; cpt++)
      shape_quad(loc_cpts[cpt], &shape_cpts(cpt,0), nNodes);

    for (uint spt = 0; spt < nSpts; spt++)
      shape_quad(loc_spts[spt], &shape_spts(spt,0), nNodes);

    // Setting unit normal vector in the parent domain
    for (uint fpt = 0; fpt < nFpts; fpt++)
    {
      shape_quad(loc_fpts[fpt], &shape_fpts(fpt,0), nNodes);

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
    loc_ppts = getLocPpts(HEX,order,params->sptsTypeQuad);
    loc_cpts = getLocSpts(HEX,order,string("Equidistant"));

    dshape_hex(loc_spts, dshape_spts, nNodes);
    dshape_hex(loc_fpts, dshape_fpts, nNodes);
    dshape_hex(loc_cpts, dshape_cpts, nNodes);

    for (uint ppt = 0; ppt < nPpts; ppt++)
      shape_hex(loc_ppts[ppt], &shape_ppts(ppt,0), nNodes);

    for (uint cpt = 0; cpt < nCpts; cpt++)
      shape_hex(loc_cpts[cpt], &shape_cpts(cpt,0), nNodes);

    for (uint spt = 0; spt < nSpts; spt++)
      shape_hex(loc_spts[spt], &shape_spts(spt,0), nNodes);

    for (uint fpt = 0; fpt < nFpts; fpt++)
    {
      shape_hex(loc_fpts[fpt], &shape_fpts(fpt,0), nNodes);

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
  /* Intermediate residuals for Runge-Kutta time integration */

  for (int step=0; step<nRKSteps-1; step++) {
    params->rkTime = params->time + params->RKa[step]*params->dt;

    if (step == 0 && params->dtType != 0) calcDt();

    moveMesh(step);

    if (step == 0) copyUspts_U0(); // Store starting values for RK method

    calcResidual(step);

    timeStepA(step, params->RKa[step+1], PMG_Source);
  }

  /* Final Runge-Kutta time advancement step */

  params->rkTime = params->time + params->RKa[nRKSteps-1]*params->dt;

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  if (params->timeType < 5)
  {
    // 'Normal' RK time-stepping: Reset solution to initial-stage values
    if (nRKSteps>1)
      copyU0_Uspts();
    else if (params->dtType != 0)
      calcDt();

    for (int step=0; step<nRKSteps; step++)
      timeStepB(step, PMG_Source);
  }
  else
  {
    // Jameson-style RK update
    timeStepA(nRKSteps-1, params->RKa[nRKSteps], PMG_Source);
  }

  params->time += params->dt;
}


void solver::calcResidual(int step)
{
  if (nEles == 0) return;

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
  if (params->iter == params->initIter+1) {
#pragma omp parallel for
    for (uint i = 0; i < nEles; i++)
      eles[i]->calcWaveSpFpts();
  }

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

void solver::timeStepA(int step, double RKval, bool PMG_Source)
{
  if (PMG_Source)
  {
    /* --- Include PMG Source Term --- */
#pragma omp parallel for collapse(3)
    for (uint spt = 0; spt < nSpts; spt++) {
      for (uint e = 0; e < nEles; e++) {
        for (uint k = 0; k < nFields; k++) {
          if (params->dtType != 2)
            U_spts(spt,e,k) = U0(spt,e,k) - RKval * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt,e,k) = U0(spt,e,k) - RKval * (divF_spts[step](spt,e,k) + src_spts(spt,e,k)) / detJac_spts(spt,e) * eles[e]->dt;
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
            U_spts(spt, e, k) = U0(spt,e,k) - RKval * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * params->dt;
          else
            U_spts(spt, e, k) = U0(spt,e,k) - RKval * divF_spts[step](spt,e,k) / detJac_spts(spt,e) * eles[e]->dt;
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
  //! TODO: Re-implement
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


void solver::extrapolateUPpts(void)
{
  int m = nPpts;
  int n = nEles * nFields;
  int k = nSpts;

#pragma omp parallel for collapse(2)
  for (int spt = 0; spt < nSpts; spt++) {
    for (int ele = 0; ele < nEles; ele++) {
      if (params->equation == ADVECTION_DIFFUSION) {
        V_spts(spt,ele,0) = U_spts(spt,ele,0);
      }
      else if (params->equation == NAVIER_STOKES) {
        double rho = U_spts(spt,ele,0);
        double u = U_spts(spt,ele,1) / rho;
        double v = U_spts(spt,ele,2) / rho;
        double w = 0;
        double vMagSq = u*u + v*v;
        if (nDims == 3) {
          w = U_spts(spt,ele,3) / rho;
          vMagSq += w*w;
          V_spts(spt,ele,3) = w;
        }
        V_spts(spt,ele,0) = rho;
        V_spts(spt,ele,1) = u;
        V_spts(spt,ele,2) = v;
        V_spts(spt,ele,nDims+1) = (params->gamma-1)*(U_spts(spt,ele,nDims+1) - 0.5*rho*vMagSq);
      }
    }
  }

  auto &A = opers[order].opp_spts_to_ppts(0,0);
  auto &B = V_spts(0,0,0);
  auto &C = V_ppts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif
}

void solver::extrapolateGridVelPpts(void)
{
  int m = nPpts;
  int n = nEles * nDims;
  int k = nNodes;

  auto &A = shape_ppts(0,0);
  auto &B = gridV_mpts(0,0,0);
  auto &C = gridV_ppts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A, k, &B, n, 0.0, &C, n);
#endif
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
  double tempF[3][5];
#pragma omp parallel for collapse(2) private(tempF)
  for (uint spt = 0; spt < nSpts; spt++) {
    for (uint e = 0; e < nEles; e++) {
      inviscidFlux(&U_spts(spt,e,0), tempF, params);

      if (params->motion || params->viscous)
      {
        /* --- Transformed later - just copy over --- */
        for (uint dim = 0; dim < nDims; dim++)
          for (uint k = 0; k < nFields; k++)
            F_spts(dim,spt,e,k) = tempF[dim][k];
      }
      else
      {
        /* --- Transform back to reference domain --- */
        for (uint dim1 = 0; dim1 < nDims; dim1++) {
          for (uint k = 0; k < nFields; k++) {
            F_spts(dim1,spt,e,k) = 0.;
            for (uint dim2 = 0; dim2 < nDims; dim2++) {
              F_spts(dim1,spt,e,k) += JGinv_spts(spt,e,dim1,dim2)*tempF[dim2][k];
            }
          }
        }
      }
    }
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
  double tempF[3][5];
#pragma omp parallel for collapse(2)
  for (uint spt = 0; spt < nSpts; spt++) {
    for (uint e = 0; e < nEles; e++) {
      for (uint dim = 0; dim < nDims; dim++)
        for (uint k = 0; k < nFields; k++)
          tempDU(dim,k) = dU_spts(dim,spt,e,k);

      if (params->equation == NAVIER_STOKES)
        viscousFlux(&U_spts(spt,e,0), tempDU, tempF, params);
      else
        viscousFluxAD(tempDU, tempF, params);

      /* Add physical inviscid flux at spts */
      for (uint dim = 0; dim < nDims; dim++)
        for (uint k = 0; k < nFields; k++)
          tempF[dim][k] += F_spts(dim,spt,e,k);

      /* --- Transform back to reference domain --- */
      for (uint dim1 = 0; dim1 < nDims; dim1++)
        for (uint k = 0; k < nFields; k++)
          F_spts(dim1,spt,e,k) = 0.;

      for (uint dim1 = 0; dim1 < nDims; dim1++) {
        for (uint k = 0; k < nFields; k++) {
          for (uint dim2 = 0; dim2 < nDims; dim2++) {
            F_spts(dim1,spt,e,k) += JGinv_spts(spt,e,dim1,dim2)*tempF[dim2][k];
          }
        }
      }
    }
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
#pragma omp parallel for collapse(2)
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
#pragma omp parallel for collapse(2) private(Jacobian)
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

  auto &A0 = opers[order].opp_grad_spts[0](0,0);
  auto &B0 = F_spts(0,0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A0, k, &B0, n, 0.0, &C, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
              1.0, &A0, k, &B0, n, 0.0, &C, n);
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
  int n = nEles * nFields;
  int k = nFpts;

  auto &B = dUc_fpts(0,0,0);

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
#pragma omp parallel for collapse(2)
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
#pragma omp parallel for collapse(2)
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
      Geo->processBlanks(eles,faces,mpiFaces,overFaces,this);
      Geo->processUnblanks(eles,faces,mpiFaces,overFaces,this);
      if (params->projection)
      {
        /* --- Use LGP with supermeshing; NOTE: For linear shape funcs only --- */
        OComm->matchUnblankCells(eles,Geo->unblankCells,Geo->eleMap,params->quadOrder);
        OComm->performGalerkinProjection(eles,opers,Geo->eleMap,order);
      }
      else
      {
        /* --- Use collocation projection; Use for nonlinear shape funcs --- */
        OComm->setupFringeCellPoints(eles,Geo->unblankCells,Geo->eleMap);
        OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
        OComm->exchangeOversetData(eles,opers,Geo->eleMap);
        OComm->transferEleData(eles,Geo->unblankCells,Geo->eleMap);
      }
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

    if (gridID == 0 || Geo->unblankCells.size()>0) {
      /*! WARNING: ONLY FOR CERTAIN MOTION TYPES! !*/
      updatePosSptsFpts();

      updateGridVSptsFpts();

      updateTransforms();
    }
    Geo->updateADT();

    if (params->oversetMethod != 2) {
      if (params->nDims==2)
        getBoundingBox(Geo->xv,Geo->minPt,Geo->maxPt);
      int nPtsFace = order+1;
      if (nDims == 3) nPtsFace *= order+1;
      OComm->setupOverFacePoints(overFaces,nPtsFace);
      OComm->matchOversetPoints(eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
    }

    if (step==0 && Geo->unblankCells.size() > 0) {
      extrapolateU();
      for (auto &e:eles)
        e->calcWaveSpFpts();
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
  if (nEles == 0) return;

  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < Geo->nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        if (Geo->eleMap[e] >= 0)
          nodes(npt, Geo->eleMap[e], dim) = Geo->xv[Geo->c2v(e,npt)][dim];

  int ms = nSpts;
  int mf = nFpts;
  int mp = nPpts;
  int mc = nCpts;
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

  auto &Ap = shape_ppts(0,0);
  auto &Cp = pos_ppts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mp, n, k,
              1.0, &Ap, k, &B, n, 0.0, &Cp, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mp, n, k,
              1.0, &Ap, k, &B, n, 0.0, &Cp, n);
#endif

  auto &Ac = shape_cpts(0,0);
  auto &Cc = pos_cpts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mp, n, k,
              1.0, &Ap, k, &B, n, 0.0, &Cp, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mc, n, k,
              1.0, &Ac, k, &B, n, 0.0, &Cc, n);
#endif

  /* Initialize storage of moving node positions */
  if (params->motion)
  {
    nodesRK = nodes;
  }
}

void solver::updatePosSptsFpts(void)
{
  if (nEles == 0) return;

#pragma omp parallel for collapse(3)
  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < Geo->nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        if (Geo->eleMap[e] >= 0)
          nodesRK(npt, Geo->eleMap[e], dim) = Geo->xv[Geo->c2v(e,npt)][dim];

  int ms = nSpts;
  int mf = nFpts;
  int mp = nPpts;
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

  auto &Ap = shape_ppts(0,0);
  auto &Cp = pos_ppts(0,0,0);
#ifdef _OMP
  omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mp, n, k,
              1.0, &Ap, k, &B, n, 0.0, &Cp, n);
#else
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mp, n, k,
              1.0, &Ap, k, &B, n, 0.0, &Cp, n);
#endif
}

void solver::updateGridVSptsFpts(void)
{
  if (nEles == 0) return;

#pragma omp parallel for collapse(3)
  for (uint npt = 0; npt < nNodes; npt++)
    for (uint e = 0; e < Geo->nEles; e++)
      for (uint dim = 0; dim < nDims; dim++)
        if (Geo->eleMap[e] >= 0)
          gridV_mpts(npt, Geo->eleMap[e], dim) = Geo->gridVel(Geo->c2v(e,npt),dim);

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

void solver::calcCSCMetrics(void)
{
  /* --- Step 1: Calulcate Inner Terms --- */

  Array<double,4> gradPos(nDims, nCpts, nEles, nDims);

  int mc = nCpts;
  int kc = nNodes;
  int nc = nEles * nDims;

  for (int dim = 0; dim < nDims; dim++) {
    auto &A = dshape_cpts(dim, 0, 0);
    auto &B = nodes(0, 0, 0);
    auto &C = gradPos(dim, 0, 0, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mc, nc, kc,
                1.0, &A, kc, &B, nc, 0.0, &C, nc);
  }

  if (nDims == 3)
  {
    // Handy cross-product indices
    int cross1[3] = {1,2,0};
    int cross2[3] = {2,0,1};

    //              xi,eta,zeta             x,y,z
    Array<double,4> TA(nDims, nCpts, nEles, nDims); // TA for xi_x, xi_y, xi_z, eta_x, eta_y, ...
    Array<double,4> TB(nDims, nCpts, nEles, nDims);
    for (int d1 = 0; d1 < 3; d1++) {
      for (int cpt = 0; cpt < nCpts; cpt++) {
        for (int e = 0; e < nEles; e++) {
          for (int d2 = 0; d2 < 3; d2++) {
            TA(d1, cpt, e, d2) = gradPos(cross1[d1],cpt,e,cross1[d2]) * pos_cpts(cpt,e,cross2[d2]) - gradPos(cross1[d1],cpt,e,cross2[d2]) * pos_cpts(cpt,e,cross1[d2]);
            TB(d1, cpt, e, d2) = gradPos(cross2[d1],cpt,e,cross1[d2]) * pos_cpts(cpt,e,cross2[d2]) - gradPos(cross2[d1],cpt,e,cross2[d2]) * pos_cpts(cpt,e,cross1[d2]);
          }
        }
      }
    }

    /* --- Step 2.1: Calulcate Jacobian Entries at spts --- */

    int ms = nSpts;
    int ks = nCpts;
    int ns = nEles * nDims;
    Array<double,3> gradTA(nSpts, nEles, nDims);
    Array<double,3> gradTB(nSpts, nEles, nDims);
    for (int d1 = 0; d1 < 3; d1++) {
      auto &AA = opers[order].gradCpts_spts[cross2[d1]](0, 0);
      auto &BA = TA(d1, 0, 0, 0);
      auto &CA = gradTA(0, 0, 0);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, ns, ks,
                  0.5, &AA, ks, &BA, ns, 0.0, &CA, ns);

      auto &AB = opers[order].gradCpts_spts[cross1[d1]](0, 0);
      auto &BB = TB(d1, 0, 0, 0);
      auto &CB = gradTB(0, 0, 0);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, ns, ks,
                  0.5, &AB, ks, &BB, ns, 0.0, &CB, ns);

      /*! !! ALTERNATIVE - REQUIRES RE-ARRANGING JGINV_SPTS !!
     * Structure as JGinv_spts(d1,spt,e,d2)
     * do dgemm on JGinv_spts(d1,0,0,0) and set beta = -1.0 for second call
     auto &C = JGinv_spts(d1,0,0,0);
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, ns, ks,
                 1.0, &AA, ks, &BA, ns, 0.0, &C, ns);
     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ms, ns, ks,
                 1.0, &AB, ks, &BB, ns, -1.0, &C, ns);
     */

      for (int spt = 0; spt < nSpts; spt++) {
        for (int e = 0; e < nEles; e++) {
          for (int d2 = 0; d2 < 3; d2++) {
            JGinv_spts(spt,e,d1,d2) = gradTA(spt,e,d2) - gradTB(spt,e,d2);
          }
        }
      }
    }

    /* --- Step 2.2: Calulcate Jacobian Entries at fpts --- */

    int mf = nFpts;
    int kf = nCpts;
    int nf = nEles * nDims;
    gradTA.setup(nFpts, nEles, nDims);
    gradTB.setup(nFpts, nEles, nDims);
    for (int d1 = 0; d1 < 3; d1++) {
      auto &AA = opers[order].gradCpts_fpts[cross2[d1]](0, 0);
      auto &BA = TA(d1, 0, 0, 0);
      auto &CA = gradTA(0, 0, 0);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, nf, kf,
                  1.0, &AA, kf, &BA, nf, 0.0, &CA, nf);

      auto &AB = opers[order].gradCpts_fpts[cross1[d1]](0, 0);
      auto &BB = TB(d1, 0, 0, 0);
      auto &CB = gradTB(0, 0, 0);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mf, nf, kf,
                  1.0, &AB, kf, &BB, nf, 0.0, &CB, nf);

      for (int fpt = 0; fpt < nFpts; fpt++) {
        for (int e = 0; e < nEles; e++) {
          for (int d2 = 0; d2 < 3; d2++) {
            JGinv_fpts(fpt,e,d1,d2) = gradTA(fpt,e,d2) - gradTB(fpt,e,d2);
          }
        }
      }
    }
  }

  /* --- Step 3: Calculate detJac and normals --- */

  for (int spt = 0; spt < nSpts; spt++) {
    for (int e = 0; e < nEles; e++) {
      double rx = JGinv_spts(spt,e,0,0);   double ry = JGinv_spts(spt,e,0,1);   double rz = JGinv_spts(spt,e,0,2);
      double sx = JGinv_spts(spt,e,1,0);   double sy = JGinv_spts(spt,e,1,1);   double sz = JGinv_spts(spt,e,1,2);
      double tx = JGinv_spts(spt,e,2,0);   double ty = JGinv_spts(spt,e,2,1);   double tz = JGinv_spts(spt,e,2,2);
      detJac_spts(spt,e) = sqrt((rx*(sy*tz - sz*ty) - ry*(sx*tz - sz*tx) + rz*(sx*ty - sy*tx))); // |JGinv| = |J|^(nDims-1)
      /// DEBUGGING CSC METRICS
      cout.precision(6);
      cout << fixed;
//      cout << " spt " << spt << ", e " << e << ", dj = " << detJac_spts(spt,e) << endl;
      if (e == 3) {
      cout << "csc rx: " << setw(10) << rx << ", ry: " << setw(10) << ry << ", rz: " << setw(10) << rz << endl;
      cout << "csc sx: " << setw(10) << sx << ", sy: " << setw(10) << sy << ", sz: " << setw(10) << sz << endl;
      cout << "csc tx: " << setw(10) << tx << ", ty: " << setw(10) << ty << ", tz: " << setw(10) << tz << endl;
      }
    }
  }

  for (int fpt = 0; fpt < nFpts; fpt++) {
    for (int e = 0; e < nEles; e++) {
      double rx = JGinv_fpts(fpt,e,0,0);   double ry = JGinv_fpts(fpt,e,0,1);   double rz = JGinv_fpts(fpt,e,0,2);
      double sx = JGinv_fpts(fpt,e,1,0);   double sy = JGinv_fpts(fpt,e,1,1);   double sz = JGinv_fpts(fpt,e,1,2);
      double tx = JGinv_fpts(fpt,e,2,0);   double ty = JGinv_fpts(fpt,e,2,1);   double tz = JGinv_fpts(fpt,e,2,2);
      detJac_fpts(fpt,e) = sqrt(rx*(sy*tz - sz*ty) - ry*(sx*tz - sz*tx) + rz*(sx*ty - sy*tx)); // |JGinv| = |J|^(nDims-1)

      /* --- Calculate outward unit normal vector at flux point --- */
      // Transform face normal from reference to physical space [JGinv .dot. tNorm]
      for (uint dim1 = 0; dim1 < nDims; dim1++) {
        norm_fpts(fpt,e,dim1) = 0.;
        for (uint dim2 = 0; dim2 < nDims; dim2++) {
          norm_fpts(fpt,e,dim1) += JGinv_fpts(fpt,e,dim2,dim1) * tNorm_fpts(fpt,dim2);
        }
      }

      // Store magnitude of face normal (equivalent to face area in finite-volume land)
      dA_fpts(fpt,e) = 0;
      for (uint dim = 0; dim < nDims; dim++)
        dA_fpts(fpt,e) += norm_fpts(fpt,e,dim)*norm_fpts(fpt,e,dim);
      dA_fpts(fpt,e) = sqrt(dA_fpts(fpt,e));

      // Normalize
      if (std::fabs(dA_fpts(fpt,e)) < 1e-10) {
        dA_fpts(fpt,e) = 0.;
        for (uint dim = 0; dim < nDims; dim++)
          norm_fpts(fpt,e,dim) = 0;
      }
      else {
        for (uint dim = 0; dim < nDims; dim++)
          norm_fpts(fpt,e,dim) /= dA_fpts(fpt,e);
      }
    }
  }
}

void solver::calcTransforms(void)
{
  /* --- Calculate Transformation at Solution Points --- */
#pragma omp parallel for collapse(2)
  for (uint spt = 0; spt < nSpts; spt++)
  {
    for (uint e = 0; e < nEles; e++)
    {
      for (uint dim1 = 0; dim1 < nDims; dim1++)
        for (uint dim2 = 0; dim2 < nDims; dim2++)
          Jac_spts(spt,e,dim1,dim2) = 0.;

      if (params->motion != 0)
      {
        for (uint i = 0; i < nNodes; i++)
          for (uint dim1 = 0; dim1 < nDims; dim1++)
            for (uint dim2 = 0; dim2 < nDims; dim2++)
              Jac_spts(spt,e,dim1,dim2) += dshape_spts(dim2,spt,i)*nodesRK(i,e,dim1);
      }
      else
      {
        for (uint i = 0; i < nNodes; i++)
          for (uint dim1 = 0; dim1 < nDims; dim1++)
            for (uint dim2 = 0; dim2 < nDims; dim2++)
              Jac_spts(spt,e,dim1,dim2) += dshape_spts(dim2,spt,i)*nodes(i,e,dim1);
      }

      if (nDims == 2) {
        // Determinant of transformation matrix
        detJac_spts(spt,e) = Jac_spts(spt,e,0,0)*Jac_spts(spt,e,1,1)-Jac_spts(spt,e,1,0)*Jac_spts(spt,e,0,1);
        // Inverse of transformation matrix (times its determinant)
        JGinv_spts(spt,e,0,0) = Jac_spts(spt,e,1,1);  JGinv_spts(spt,e,0,1) =-Jac_spts(spt,e,0,1);
        JGinv_spts(spt,e,1,0) =-Jac_spts(spt,e,1,0);  JGinv_spts(spt,e,1,1) = Jac_spts(spt,e,0,0);
      }
      else if (nDims == 3) {
        double xr = Jac_spts(spt,e,0,0);   double xs = Jac_spts(spt,e,0,1);   double xt = Jac_spts(spt,e,0,2);
        double yr = Jac_spts(spt,e,1,0);   double ys = Jac_spts(spt,e,1,1);   double yt = Jac_spts(spt,e,1,2);
        double zr = Jac_spts(spt,e,2,0);   double zs = Jac_spts(spt,e,2,1);   double zt = Jac_spts(spt,e,2,2);
        detJac_spts(spt,e) = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);

        JGinv_spts(spt,e,0,0) = ys*zt - yt*zs;  JGinv_spts(spt,e,0,1) = xt*zs - xs*zt;  JGinv_spts(spt,e,0,2) = xs*yt - xt*ys;
        JGinv_spts(spt,e,1,0) = yt*zr - yr*zt;  JGinv_spts(spt,e,1,1) = xr*zt - xt*zr;  JGinv_spts(spt,e,1,2) = xt*yr - xr*yt;
        JGinv_spts(spt,e,2,0) = yr*zs - ys*zr;  JGinv_spts(spt,e,2,1) = xs*zr - xr*zs;  JGinv_spts(spt,e,2,2) = xr*ys - xs*yr;
      }
      if (detJac_spts(spt,e)<0) FatalError("Negative Jacobian at solution points.");
    }
  }

  /* --- Calculate Transformation at Flux Points --- */
#pragma omp parallel for collapse(2)
  for (uint fpt = 0; fpt < nFpts; fpt++)
  {
    for (uint e = 0; e < nEles; e++)
    {
      for (uint dim1 = 0; dim1 < nDims; dim1++)
        for (uint dim2 = 0; dim2 < nDims; dim2++)
          Jac_fpts(fpt,e,dim1,dim2) = 0.;

      // Calculate transformation Jacobian matrix - [dx/dr, dx/ds; dy/dr, dy/ds]
      if (params->motion != 0)
      {
        for (uint i = 0; i < nNodes; i++)
          for (uint dim1 = 0; dim1 < nDims; dim1++)
            for (uint dim2 = 0; dim2 < nDims; dim2++)
              Jac_fpts(fpt,e,dim1,dim2) += dshape_fpts(dim2,fpt,i)*nodesRK(i,e,dim1);
      }
      else
      {
        for (uint i = 0; i < nNodes; i++)
          for (uint dim1 = 0; dim1 < nDims; dim1++)
            for (uint dim2 = 0; dim2 < nDims; dim2++)
              Jac_fpts(fpt,e,dim1,dim2) += dshape_fpts(dim2,fpt,i)*nodes(i,e,dim1);
      }

      if (nDims == 2) {
        detJac_fpts(fpt,e) = Jac_fpts(fpt,e,0,0)*Jac_fpts(fpt,e,1,1)-Jac_fpts(fpt,e,1,0)*Jac_fpts(fpt,e,0,1);
        // Inverse of transformation matrix (times its determinant)
        JGinv_fpts(fpt,e,0,0) = Jac_fpts(fpt,e,1,1);  JGinv_fpts(fpt,e,0,1) =-Jac_fpts(fpt,e,0,1);
        JGinv_fpts(fpt,e,1,0) =-Jac_fpts(fpt,e,1,0);  JGinv_fpts(fpt,e,1,1) = Jac_fpts(fpt,e,0,0);
      }
      else if (nDims == 3) {
        double xr = Jac_fpts(fpt,e,0,0);   double xs = Jac_fpts(fpt,e,0,1);   double xt = Jac_fpts(fpt,e,0,2);
        double yr = Jac_fpts(fpt,e,1,0);   double ys = Jac_fpts(fpt,e,1,1);   double yt = Jac_fpts(fpt,e,1,2);
        double zr = Jac_fpts(fpt,e,2,0);   double zs = Jac_fpts(fpt,e,2,1);   double zt = Jac_fpts(fpt,e,2,2);
        detJac_fpts(fpt,e) = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);
        // Inverse of transformation matrix (times its determinant)
        JGinv_fpts(fpt,e,0,0) = ys*zt - yt*zs;  JGinv_fpts(fpt,e,0,1) = xt*zs - xs*zt;  JGinv_fpts(fpt,e,0,2) = xs*yt - xt*ys;
        JGinv_fpts(fpt,e,1,0) = yt*zr - yr*zt;  JGinv_fpts(fpt,e,1,1) = xr*zt - xt*zr;  JGinv_fpts(fpt,e,1,2) = xt*yr - xr*yt;
        JGinv_fpts(fpt,e,2,0) = yr*zs - ys*zr;  JGinv_fpts(fpt,e,2,1) = xs*zr - xr*zs;  JGinv_fpts(fpt,e,2,2) = xr*ys - xs*yr;
      }

      /* --- Calculate outward unit normal vector at flux point --- */
      // Transform face normal from reference to physical space [JGinv .dot. tNorm]
      for (uint dim1 = 0; dim1 < nDims; dim1++) {
        norm_fpts(fpt,e,dim1) = 0.;
        for (uint dim2 = 0; dim2 < nDims; dim2++) {
          norm_fpts(fpt,e,dim1) += JGinv_fpts(fpt,e,dim2,dim1) * tNorm_fpts(fpt,dim2);
        }
      }

      // Store magnitude of face normal (equivalent to face area in finite-volume land)
      dA_fpts(fpt,e) = 0;
      for (uint dim = 0; dim < nDims; dim++)
        dA_fpts(fpt,e) += norm_fpts(fpt,e,dim)*norm_fpts(fpt,e,dim);
      dA_fpts(fpt,e) = sqrt(dA_fpts(fpt,e));

      // Normalize
      // If we have a collapsed edge, the dA will be 0, so just set the normal to 0
      // (A normal vector at a point doesn't make sense anyways)
      if (std::fabs(dA_fpts(fpt,e)) < 1e-10) {
        dA_fpts(fpt,e) = 0.;
        for (uint dim = 0; dim < nDims; dim++)
          norm_fpts(fpt,e,dim) = 0;
      }
      else {
        for (uint dim = 0; dim < nDims; dim++)
          norm_fpts(fpt,e,dim) /= dA_fpts(fpt,e);
      }
    }
  }
}

void solver::updateTransforms(void)
{
  if (params->motion != 4)
  {
    calcTransforms();
  }
  else
  {
    for (auto &ic:Geo->unblankCells)
      eles[Geo->eleMap[ic]]->calcTransforms(true);
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

    Geo->processBlanks(eles,faces,mpiFaces,overFaces,this);
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
      int nPtsFace = order+1;
      if (nDims==3) nPtsFace *= order+1;
      if (params->oversetMethod != 2)
        OComm->setupOverFacePoints(overFaces, nPtsFace);
      else
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
  vector<double> LpErr(params->nFields);

  if (params->errorNorm == 0)
  {
    /* Integrate error and solution over each element */

    auto wts = getQptWeights(order, nDims);
    vector<double> intU(params->nFields);
    for (uint ic = 0; ic < eles.size(); ic++) {
      auto tmpUE = eles[ic]->calcError();
      for (uint spt = 0; spt < nSpts; spt++) {
        for (uint k = 0; k < nFields; k++) {
          LpErr[k] += tmpUE(spt,k) * wts[spt] * detJac_spts(spt, ic);
          intU[k] += U_spts(spt,ic,k) * wts[spt] * detJac_spts(spt, ic);
        }
      }
    }

    for (uint k = 0; k < nFields; k++)
      LpErr[k] = abs(LpErr[k] - intU[k]);
  }
  else
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

  #ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
  #else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B, n, 0.0, &C, n);
  #endif

    n = nEles;
    auto &B1 = detJac_spts(0,0);
    auto &C1 = detJac_qpts(0,0);
  #ifdef _OMP
    omp_blocked_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B1, n, 0.0, &C1, n);
  #else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                1.0, &A, k, &B1, n, 0.0, &C1, n);
  #endif

    /* Integrate error over each element */

    for (uint ic = 0; ic < eles.size(); ic++) {
      for (uint qpt = 0; qpt < nQpts; qpt++) {
        auto tmpErr = calcError(&U_qpts(qpt, ic, 0), eles[ic]->calcPos(qpts[qpt]), params);
        for (uint j = 0; j < params->nFields; j++)
          LpErr[j] += tmpErr[j] * wts[qpt] * detJac_qpts(qpt, ic);
      }
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
  //! TODO: Re-implement
//#pragma omp parallel for
//  for (uint i=0; i<eles.size(); i++) {
//    eles[i]->sensor = opers[order].shockCaptureInEle(eles[i]->U_spts,params->threshold);
//  }
}
