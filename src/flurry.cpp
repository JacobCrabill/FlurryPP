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

  int iter;

  cout << "< insert ASCII art here >" << endl;

  if (argc<2) FatalError("No input file specified.");

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  /* Setup the mesh and connectivity for the simulation */
  Geo.setup(&params);

  /* Setup the solver, all elements and faces, and all FR operators for computation */
  Solver.setup(&params,&Geo);

  /* Apply the initial condition */
  Solver.initializeSolution();

  /* Write initial data file */
  writeData(&Solver,&params,0);

  /* --- Calculation Loop --- */
  for (iter=params.initIter; iter<params.iterMax; iter++) {

    Solver.calcResidual();

    if ((iter+1)%params.monitor_res_freq == 0) writeResidual(&Solver,&params,iter);
    if ((iter+1)%params.plot_freq == 0) writeData(&Solver,&params,iter);

  }
}
