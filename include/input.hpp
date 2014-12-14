/*!
 * \file input.hpp
 * \brief Header file for input reader
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Fux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */
#include "../include/global.hpp"

#include <string>
#include <fstream>

class input
{
public:
  //! Default initialization
  input();

  void read_input_file(char *filename);

  /* --- Member Variables --- */
  int equation;  //! {0 | Advection/diffusion} {1 | Euler/Navier-Stokes}
  int viscous;   //! {0 | No diffusion/Viscosity} {1 | Diffusion/Viscosity}
  int order;
  int riemann_type;
  int ic_type;
  int test_case;

  /* --- Simulation Run Parameters --- */
  double dt;
  int iterMax;
  int plot_freq;
  int restart_freq;


  /* --- Boundary Condition Parameters --- */
  double Uinf;
  double Vinf;
  double rhoinf;
  double Pinf;

  string spt_type_tri;  //! Legendre, Lobatto, ...
  string spt_type_quad;

private:

};
