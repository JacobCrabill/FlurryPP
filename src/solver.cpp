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

void solver::initialize(input *params)
{
  this->params = params;
}

void solver::calcResidual(void)
{
  extrapolateU();

  calcInviscidFlux_spts();

  calcInviscidFlux_faces();

  if (params.viscous) {
    calcGradU_spts();

    correctU();

    extrapolateGradU();

    calcViscousFlux_faces();
  }

  if (params.motion) {
    // Use non-conservative chain-rule formulation (Liang-Miyaji)
    calcGradF_spts();
  }else{
    // Standard conservative form
    calcDivF_spts();
  }

  correctFlux();
}

void solver::extrapolateU(void)
{
  for (auto& e: eles) {
    opers[e.eType][e.order].apply_spts_fpts(e.U_spts,e.U_fpts);
  }
}

void solver::calcInviscidFlux_spts(void)
{
  for (auto& e:eles) {
    e.calcInviscidFlux();
  }
}

void solver::calcViscousFlux_spts(void)
{
  for (auto& e:eles) {
    e.calcInviscidFlux();
  }
}

void solver::calc_grad_spts(void)
{
  for (auto& e: eles) {
    opers[e.eType][e.order].apply_grad_spts(e.U_spts,e.dU_spts);
  }
}

void solver::apply_oper(matrix<double> &op, matrix<double> & A, matrix<double> &B)
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
    for (auto& p: polyOrders[eTypes]) {
      opers[e][p].setupOperators(e,p);
    }
  }
}
