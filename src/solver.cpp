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

solver::solver()
{
}

void solver::setup(input *params, geo *Geo)
{
  this->params = params;
  this->Geo = Geo;

  /* Setup the FR elements & faces which will be computed on */
  Geo->setupElesFaces(eles,faces,bounds);

  /* Setup the FR operators for computation */
  setupOperators();
}

void solver::calcResidual(void)
{
  extrapolateU();

  calcInviscidFlux_spts();

  extrapolateNormalFlux();

  /* Inviscid Common Flux */
  calcInviscidFlux_faces();

  calcInviscidFlux_bounds();

  if (params->viscous) {
    calcGradU_spts();

    correctU();

    extrapolateGradU();

    /* Viscous Common Flux */
    calcViscousFlux_faces();

    calcViscousFlux_bounds();
  }

  if (params->motion) {
    /* Use non-conservative chain-rule formulation (Liang-Miyaji) */
    calcGradF_spts();
  }else{
    /* Standard conservative form */
    calcDivF_spts();
  }

  /* Extrapolate total flux to flux points & dot with normal */
  //extrapolateNormalFlux();

  correctDivFlux();
}

void solver::timeStep(void)
{
  for (auto& e:eles) {
    e.timeStep();
  }
}

void solver::extrapolateU(void)
{
  for (auto& e:eles) {
    opers[e.eType][e.order].applySptsFpts(e.U_spts,e.U_fpts);
  }
}

void solver::calcInviscidFlux_spts(void)
{
  for (auto& e:eles) {
    e.calcInviscidFlux_spts();
  }
}

void solver::calcInviscidFlux_faces()
{
  for (auto& F:faces) {
    F.calcInviscidFlux();
  }
}

void solver::calcInviscidFlux_bounds()
{
  for (auto& B:bounds) {
    B.calcInviscidFlux();
  }
}

void solver::calcViscousFlux_spts(void)
{
  for (auto& e:eles) {
    e.calcInviscidFlux_spts();
  }
}

void solver::calcViscousFlux_faces()
{

}

void solver::calcViscousFlux_bounds()
{

}

void solver::calcGradF_spts()
{
  for (auto& e: eles) {
    opers[e.eType][e.order].applyGradFSpts(e.F_spts,e.dF_spts);
  }
}

void solver::calcDivF_spts()
{
  for (auto& e: eles) {
    opers[e.eType][e.order].applyDivFSpts(e.F_spts,e.divF_spts);
  }
}

void solver::extrapolateNormalFlux(void)
{
  for (auto& e: eles) {
    opers[e.eType][e.order].applyExtrapolateFn(e.F_spts,e.tNorm_fpts,e.Fn_fpts);
  }
}

void solver::correctDivFlux()
{
  for (auto& e: eles) {
    opers[e.eType][e.order].applyCorrectDivF(e.dFn_fpts,e.divF_spts);
  }

}

void solver::calcGradU_spts(void)
{
  for (auto& e: eles) {
    opers[e.eType][e.order].applyGradSpts(e.U_spts,e.dU_spts);
  }
}

void solver::correctU()
{

}

void solver::extrapolateGradU()
{

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
  for (auto& e:eles) {
    e.setInitialCondition();
  }
}
