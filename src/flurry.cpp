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
#include <chrono>

using namespace std::chrono;

#include "../include/flurry.hpp"

int main(int argc, char *argv[]) {
  input params;
  geo Geo;
  solver Solver;

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

  setGlobalVariables();

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  /* Setup the mesh and connectivity for the simulation */
  Geo.setup(&params);

  /* Setup the solver, all elements and faces, and all FR operators for computation */
  Solver.setup(&params,&Geo);

  /* Starting time for simulation (ignoring pre-processing) */
  high_resolution_clock::time_point initTime = high_resolution_clock::now();

  /* Apply the initial condition */
  Solver.initializeSolution();

  /* Write initial data file */
  writeData(&Solver,&params);

  /* --- Calculation Loop --- */
  for (params.iter=params.initIter+1; params.iter<=params.iterMax; params.iter++) {

    Solver.calcResidual();

    Solver.timeStep();

    if ((params.iter)%params.monitor_res_freq == 0 || params.iter==1) writeResidual(&Solver,&params);
    if ((params.iter)%params.plot_freq == 0) writeData(&Solver,&params);

  }

  // Get simulation wall time
  high_resolution_clock::time_point finalTime = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( finalTime - initTime ).count();
  double execTime = (double)duration/1000.;
  cout << setprecision(3) << "Execution time = " << execTime << "s" << endl;
}
