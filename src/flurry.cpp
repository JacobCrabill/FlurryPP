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

  clock_t initTime, finalTime;
  initTime = clock();

  cout << "  ========================================== " << endl;
  cout << "   _______   _                               " << endl;
  cout << "  |   ____| | |                              " << endl;
  cout << "  |  |___   | |  _   _   _     _     _    _  " << endl;
  cout << "  |   ___|  | | | | | | | |/| | |/| | |  | | " << endl;
  cout << "  |  |      | | | |_| | |  /  |  /  \\  \\/  / " << endl;
  cout << "  |__|      |_| \\_____/ |_|   |_|    \\    /  " << endl;
  cout << "                                      |  /   " << endl;
  cout << "                                      /_/    " << endl;
  cout << "  ----    Flux Reconstruction in C++   ----" << endl;
  cout << "  ========================================== " << endl;

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
  for (iter=params.initIter+1; iter<=params.iterMax; iter++) {

    Solver.calcResidual();

    Solver.timeStep();

    if ((iter)%params.monitor_res_freq == 0 || iter==1) writeResidual(&Solver,&params,iter);
    if ((iter)%params.plot_freq == 0) writeData(&Solver,&params,iter);

  }

  finalTime = clock()-initTime;
  printf("Execution time= %f s\n", (double) finalTime/((double) CLOCKS_PER_SEC));
}
