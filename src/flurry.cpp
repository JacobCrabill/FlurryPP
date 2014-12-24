/*!
 * \file flurry.cpp
 * \brief Main file to run a whole simulation
 *
 * Will make into a class in the future in order to interface with Python
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

#include "../include/flurry.hpp"

int main(int argc, char *argv[]) {
  input params;
  geo Geo;
  solver Solver;

  params.readInputFile(argv[1]);
  Geo.setup(&params);
  Solver.initialize(&params);
  Geo.setupElesFaces(&Solver);
}
