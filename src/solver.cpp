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
    calGradU_spts();

    correctU();

    extrapolateGradU();

    calcViscousFlux_faces();
  }
  if (params.motion) {
    // Use non-conservative chain-rule formulate (Liang-Miyaji)
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
    opers[e.eType].get_opp_spts_fpts[e.order].timesMatrix(e.U_spts,e.U_fpts);
  }
}

solver::calc_grad_spts(void)
{
  for (auto& e: eles) {
    opers[e.eType].get_opp_grad_spts[e.order].timesMatrix(e.U_spts,e.dU_spts);
  }
}

solver::apply_oper(matrix<double> &op, matrix<double> & A, matrix<double> &B)
{

}


solver::setupOperators()
{
  for (auto& e: eTypes) {
    for (auto& p: polyOrders[eTypes]) {
      //! Setup extrapolation operator
      opers[eType].setup_operators(e,p);
    }
  }
}
