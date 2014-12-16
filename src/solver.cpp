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

solver::calc_gradient()
{

}

solver::extrapolate_spts_fpts()
{
  for (auto& e: eles) {
    //! Extrapolate the solution to the flux points
    opers[e.eType].get_opp_spts_fpts[e.order].timesMatrix(e.U_spts,e.U_fpts);
  }
}

solver::apply_oper(matrix<double> &op, matrix<double> & A, matrix<double> &B)
{

}


solver::setup_operators()
{
  for (auto& e: eTypes) {
    for (auto& p: polyOrders[eTypes]) {
      //! Setup extrapolation operator
      opers[eType].setup_operators(e,p);
    }
  }
}
