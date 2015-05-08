/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include "../include/solver.hpp"

#include <omp.h>

solver::solver()
{
}

void solver::setup(input *params, geo *Geo)
{
  this->params = params;
  this->Geo = Geo;

  params->time = 0.;

  /* Setup the FR elements & faces which will be computed on */
  Geo->setupElesFaces(eles,faces,bounds);

  /* Setup the FR operators for computation */
  setupOperators();

  /* Additional Setup */

  // Time advancement setup
  switch (params->timeType) {
    case 0:
      nRKSteps = 1;
      RKb = {1};
      break;
    case 4:
      nRKSteps = 4;
      RKa = {.5, .5, 1.};
      RKb = {1./6., 1./3., 1./3., 1./6.};
      break;
    default:
      FatalError("Time-Stepping type not supported.");
  }
}

void solver::update(void)
{
  if (nRKSteps>1)
    copyUspts_U0();

  /* Intermediate residuals for Runge-Kutta time integration */

  for (int step=0; step<nRKSteps-1; step++) {

    if (step == 1)
      params->rkTime = params->time;
    else
      params->rkTime = params->time + RKa[step-1]*params->dt;

    moveMesh(step);

    calcResidual(step);

    timeStepA(step);

  }

  /* Final Runge-Kutta time advancement step */

  if (nRKSteps == 1)
    params->rkTime = params->time;
  else
    params->rkTime = params->time + params->dt;

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  if (nRKSteps>1)
    copyU0_Uspts();

  for (int step=0; step<nRKSteps; step++) {
    timeStepB(step);
  }

  params->time += params->dt;
}

void solver::calcResidual(int step)
{
  extrapolateU();

  calcInviscidFlux_spts();

  extrapolateNormalFlux();

  calcInviscidFlux_faces();

  calcInviscidFlux_bounds();

  if (params->viscous || params->motion) {

    calcGradU_spts();

  }

  if (params->viscous) {

    correctU();

    extrapolateGradU();

    calcViscousFlux_faces();

    calcViscousFlux_bounds();

  }

  calcFluxDivergence(step);

  correctDivFlux(step);
}

void solver::timeStepA(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].timeStepA(step,RKa[step]);
  }

  params->rkTime = params->time + RKa[step]*params->dt;
}

void solver::timeStepB(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].timeStepB(step,RKb[step]);
  }
}

void solver::copyUspts_U0(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].copyUspts_U0();
  }
}

void solver::copyU0_Uspts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].copyU0_Uspts();
  }
}

void solver::extrapolateU(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsFpts(eles[i].U_spts,eles[i].U_fpts);
  }
}

void solver::extrapolateUMpts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applySptsMpts(eles[i].U_spts,eles[i].U_mpts);
  }
}

void solver::calcInviscidFlux_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcInviscidFlux_spts();
  }
}

void solver::calcInviscidFlux_faces()
{
#pragma omp parallel for
  for (uint i=0; i<faces.size(); i++) {
    faces[i].calcInviscidFlux();
  }
}

void solver::calcInviscidFlux_bounds()
{
#pragma omp parallel for
  for (uint i=0; i<bounds.size(); i++) {
    bounds[i].calcInviscidFlux();
  }
}

void solver::calcViscousFlux_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].calcInviscidFlux_spts();
  }
}

void solver::calcViscousFlux_faces()
{

}

void solver::calcViscousFlux_bounds()
{

}

void solver::calcGradF_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyGradFSpts(eles[i].F_spts,eles[i].dF_spts);
  }
}

void solver::transformGradF_spts(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].transformGradF_spts(step);
    //opers[eles[i].eType][eles[i].order].applyTransformGradFSpts(eles[i].dF_spts,eles[i].JGinv_spts,eles[i].gridVel_spts);
  }
}

void solver::calcFluxDivergence(int step)
{
  if (params->motion) {

    /* Use non-conservation-form chain-rule transformation
     * (See AIAA paper 2013-0998 by Liang, Miyaji and Zhang) */

    calcGradF_spts();

    transformGradF_spts(step);

  }else{

    /* Standard conservative form */
    calcDivF_spts(step);

  }
}

void solver::calcDivF_spts(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyDivFSpts(eles[i].F_spts,eles[i].divF_spts[step]);
  }
}

void solver::extrapolateNormalFlux(void)
{
  if (params->motion)
#pragma omp parallel for
    for (uint i=0; i<eles.size(); i++) {
      opers[eles[i].eType][eles[i].order].applyExtrapolateFn(eles[i].F_spts,eles[i].norm_fpts,eles[i].Fn_fpts,eles[i].dA_fpts);
    }
  else {
    for (uint i=0; i<eles.size(); i++) {
      opers[eles[i].eType][eles[i].order].applyExtrapolateFn(eles[i].F_spts,eles[i].tNorm_fpts,eles[i].Fn_fpts,eles[i].dA_fpts);
    }
  }
}

void solver::correctDivFlux(int step)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyCorrectDivF(eles[i].dFn_fpts,eles[i].divF_spts[step]);
  }

}

void solver::calcGradU_spts(void)
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    opers[eles[i].eType][eles[i].order].applyGradSpts(eles[i].U_spts,eles[i].dU_spts);
  }
}

void solver::correctU()
{

}

void solver::extrapolateGradU()
{

}

void solver::moveMesh(int step)
{
  if (!params->motion) return;

#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].move(step);
  }
}

void solver::setupOperators()
{
  // Get all element types & olynomial orders in mesh
  for (auto& e:eles) {
    eTypes.insert(e.eType);
    polyOrders[e.eType].insert(e.order);
  }

  for (auto& e: eTypes) {
    for (auto& p: polyOrders[e]) {
      opers[e][p].setupOperators(e,p,Geo,params);
    }
  }
}

void solver::initializeSolution()
{
#pragma omp parallel for
  for (uint i=0; i<eles.size(); i++) {
    eles[i].setInitialCondition();
  }
}
