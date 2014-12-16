/*!
 * \file operators.hpp
 * \brief Header file for oper class
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
#pragma once

#include "global.hpp"
#include "mesh.hpp"

class oper
{
public:
  void setup_operators(int eleType, int order);

  //void interpolate_solution_basis(int eleType, int order, );
  //void correct_vector();
  //void correct_scalar();

  // Create a map<int,double*> (?) to get access to the correct operator
  // i.e. somthing like: div_flux_spts_tri = oper.get_oper_div[TRI]

  map<int,matrix<double>*> get_oper_div;
  map<int,matrix<double>*> get_oper_grad;
  map<int,matrix<double>*> get_oper_spts_fpts;
  map<int,matrix<double>*> get_oper_correct;

  matrix<double> opp_spts_to_fpts;
  matrix<double> opp_grad_spts;
  matrix<double> opp_correction;

private:
  mesh *Mesh;
}
