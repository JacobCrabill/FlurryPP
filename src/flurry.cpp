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
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 * Flurry++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Flurry++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Flurry++; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "flurry.hpp"

#include <valgrind/callgrind.h>

#ifndef _NO_MPI
#include <mpi.h>
#endif

#ifdef _MPI_DEBUG
#include <unistd.h>  // for getpid()
#endif

#include "funcs.hpp"
#include "multigrid.hpp"

int main(int argc, char *argv[]) {
  input params;
  solver Solver;
  multiGrid pmg;

  int rank = 0;
  int nproc = 1;
#ifndef _NO_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
  params.rank = rank;
  params.nproc = nproc;
CALLGRIND_STOP_INSTRUMENTATION;
  if (rank == 0) {
    cout << endl;
    cout << R"(  ========================================================  )" << endl;
    cout << R"(   _______   _                                     /\       )" << endl;
    cout << R"(  |   ____| | |                               __   \/   __  )" << endl;
    cout << R"(  |  |___   | |  _   _   _     _     _    _   \_\_\/\/_/_/  )" << endl;
    cout << R"(  |   ___|  | | | | | | | |/| | |/| | |  | |    _\_\/_/_    )" << endl;
    cout << R"(  |  |      | | | |_| | |  /  |  /  \  \/  /   __/_/\_\__   )" << endl;
    cout << R"(  |__|      |_| \_____/ |_|   |_|    \    /   /_/  \/  \_\  )" << endl;
    cout << R"(                                      |  /         /\       )" << endl;
    cout << R"(                                      /_/          \/       )" << endl;
    cout << R"(  ---------      Flux Reconstruction in C++      ---------  )" << endl;
    cout << R"(  ========================================================  )" << endl;
    cout << endl;
  }

  if (argc<2) FatalError("No input file specified.");

#ifdef _MPI_DEBUG
  /*// Uncomment for use with GDB or other debugger
  {
    // Useful for debugging in parallel with GDB or similar debugger
    // Sleep until a debugger is attached to rank 0 (change 'blah' once attached to continue)
    if (rank == 0) {
      int blah = 0;
      cout << "Process " << getpid() << " ready for GDB attach" << endl;
      while (blah == 0)
        sleep(5);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }*/
#endif

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  if (params.PMG)
  {
    /* Setup the P-Multigrid class if requested */
    pmg.setup(params.order,&params,Solver);
  }
  else
  {
    /* Setup the solver, grid, all elements and faces, and all FR operators for computation */
    Solver.setup(&params,params.order);

    /* Apply the initial condition */
    Solver.initializeSolution();
  }

  /* Write initial data file */
  writeData(&Solver,&params);

#ifndef _NO_MPI
  // Allow all processes to finish initial file writing before starting computation
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Start timer for simulation (ignoring pre-processing) */
  params.timer.startTimer();

  double maxTime = params.maxTime;
  int initIter = params.initIter;
  int iterMax = params.iterMax;
  int &iter = params.iter;
  iter = initIter;

  /* --- Calculation Loop --- */
  CALLGRIND_START_INSTRUMENTATION;
  while (params.iter < iterMax and params.time < maxTime) {
    iter++;

    Solver.update();

    /* If using multigrid, perform correction cycle */
    if (params.PMG)
      pmg.cycle(Solver);

    if ((iter)%params.monitorResFreq==0 or iter==initIter+1 or params.time>=maxTime) writeResidual(&Solver,&params);
    if ((iter)%params.monitorErrFreq==0 or iter==initIter+1) writeError(&Solver,&params);
    if ((iter)%params.plotFreq==0 or iter==iterMax or params.time>=maxTime) writeData(&Solver,&params);
  }
  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;

  /* Calculate the integral / L1 / L2 error for the final time */
  writeAllError(&Solver,&params);

  // Get simulation wall time
  params.timer.stopTimer();
  params.timer.showTime();

#ifndef _NO_MPI
 MPI_Finalize();
#endif

 return 0;
}
