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

#ifndef _NO_MPI
#include <mpi.h>
#endif

#ifdef _MPI_DEBUG
#include <unistd.h>  // for getpid()
#endif

#include "funcs.hpp"

int main(int argc, char *argv[]) {
  input params;
  geo Geo;
  solver Solver;

  int rank = 0;
  int nproc = 1;
#ifndef _NO_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
  params.rank = rank;
  params.nproc = nproc;

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

  setGlobalVariables();

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  bool DO_LGP_TEST = false;
  if (DO_LGP_TEST) {
    params.oversetGrids[0] = string("box1.msh");
    params.oversetGrids[1] = string("box2.msh");
  }

  /* Setup the mesh and connectivity for the simulation */
  Geo.setup(&params);

  /* Setup the solver, all elements and faces, and all FR operators for computation */
  Solver.setup(&params,&Geo);

  /* Apply the initial condition */
  Solver.initializeSolution();

#ifndef _NO_MPI
  if (DO_LGP_TEST) {
    writeData(&Solver,&params);

    vector<double> Error0, Error1, Error2, Error3;
    // Calc Error Before Interpolation:
    Error0 = Solver.integrateError();

    Geo.fringeCells.clear();
    if (Geo.gridID>0)
      for (int ic=0; ic<Geo.nEles; ic++)
        Geo.fringeCells.insert(ic);

    // Interpolate IC to 2nd grid [from box1 to box2]
    Solver.OComm->matchUnblankCells(Solver.eles,Geo.fringeCells,Geo.eleMap,10);
    Solver.OComm->performGalerkinProjection(Solver.eles,Solver.opers,Geo.eleMap);

    params.dataFileName = "LGP_Test1";
    writeData(&Solver,&params);

    // Calc Error After 1st Interpolation:
    Error1 = Solver.integrateError();

    // Setup next grid
    solver Solver2;
    geo Geo2;
    input params2 = params;
    params2.oversetGrids[0] = string("box3.msh");

    Geo2.setup(&params2);
    Solver2.setup(&params2,&Geo2);
    if (Geo.gridID>0)
      Solver2.eles = Solver.eles;

    Geo2.fringeCells.clear();
    if (Geo2.gridID==0) {
      for (int ic=0; ic<Geo2.nEles; ic++)
        Geo2.fringeCells.insert(ic);
    }

    // Interpolate IC to 2nd grid [from box2 to box3]
    Solver2.OComm->matchUnblankCells(Solver2.eles,Geo2.fringeCells,Geo2.eleMap,10);
    Solver2.OComm->performGalerkinProjection(Solver2.eles,Solver2.opers,Geo2.eleMap);

    params2.dataFileName = "LGP_Test2";
    writeData(&Solver2,&params2);

    // Calc Error After 2nd Interpolation:
    Error2 = Solver2.integrateError();


    // Setup next grid
    solver Solver3;
    geo Geo3;
    input params3 = params2;
    params3.oversetGrids[0] = string("box3.msh");
    params3.oversetGrids[1] = string("box1.msh");

    Geo3.setup(&params3);
    Solver3.setup(&params3,&Geo3);
    if (Geo.gridID==0)
      Solver3.eles = Solver2.eles;

    Geo3.fringeCells.clear();
    if (Geo3.gridID>0) {
      for (int ic=0; ic<Geo3.nEles; ic++)
        Geo3.fringeCells.insert(ic);
    }

    // Interpolate IC to 2nd grid [from box2 to box3]
    Solver3.OComm->matchUnblankCells(Solver3.eles,Geo3.fringeCells,Geo3.eleMap,10);
    Solver3.OComm->performGalerkinProjection(Solver3.eles,Solver3.opers,Geo3.eleMap);

    params3.dataFileName = "LGP_Test3";
    writeData(&Solver3,&params3);

    // Calc Error After 2nd Interpolation:
    Error3 = Solver3.integrateError();


    if (params.rank==0) cout << endl;
    double EXACT = pi*erf(5)*erf(5);
    MPI_Barrier(MPI_COMM_WORLD);

    if (Geo.gridRank==0) {
      cout << "Integral for First Grid_" << Geo.gridID << ": ";
      for (auto &val:Error0)
        cout << setprecision(16) << EXACT - val << " ";
      cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (params.rank==0) cout << endl;

    if (Geo.gridRank==0) {
      cout << "Integral for Second Grid_" << Geo.gridID << ": ";
      for (auto &val:Error1)
        cout << setprecision(16) << EXACT - val << " ";
      cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (params.rank==0) cout << endl;

    if (Geo.gridRank==0) {
      cout << "Integral for Third Grid_" << Geo.gridID << ": ";
      for (auto &val:Error2)
        cout << setprecision(16) << EXACT - val << " ";
      cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (params.rank==0) cout << endl;

    if (Geo.gridRank==0) {
      cout << "Integral for Fourth Grid_" << Geo.gridID << ": ";
      for (auto &val:Error3)
        cout << setprecision(16) << EXACT - val << " ";
      cout << endl;
    }

    MPI_Finalize();
    return 0;
  }
#endif

  /* Write initial data file */
  writeData(&Solver,&params);

#ifndef _NO_MPI
  // Allow all processes to finish initial file writing before starting computation
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Start timer for simulation (ignoring pre-processing) */
  simTimer runTime;
  runTime.startTimer();

  double maxTime = params.maxTime;
  int initIter = params.initIter;
  int iterMax = params.iterMax;
  int iter = initIter;

  /* --- Calculation Loop --- */
  while (params.iter < iterMax and params.time < maxTime) {
    iter++;

    Solver.update();

    if ((iter)%params.monitorResFreq==0 or iter==initIter+1 or params.time>=maxTime) writeResidual(&Solver,&params);
    if ((iter)%params.monitorErrFreq==0 or iter==initIter+1) writeError(&Solver,&params);
    if ((iter)%params.plotFreq==0 or iter==iterMax or params.time>=maxTime) writeData(&Solver,&params);
  }

  /* Calculate the integral / L1 / L2 error for the final time */
  writeAllError(&Solver,&params);

  // Get simulation wall time
  runTime.stopTimer();
  runTime.showTime();

  if (DO_LGP_TEST) {
    params.errorNorm = 0;
    auto err = Solver.integrateError();
    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    if (Geo.gridRank==0)
      cout << "Integrated Conservation Error for Grid " << Geo.gridID << " = " << err[0] << endl;
  } else if (params.nDims==2 && params.equation==NAVIER_STOKES) {
    params.errorNorm = 0;
    params.testCase = 0;
    vector<double> EXACT = {99.804294352079, 89.2676720089582, 44.6338359173803, 2031.69314321021}; // For L-M vortex case, 10x10 grid
    vector<double> err;
    if (params.meshType == OVERSET_MESH)
      err = Solver.integrateErrorOverset();
    else
      err = Solver.integrateError();
    cout.precision(10);
    cout.setf(ios::scientific, ios::floatfield);
    if (params.rank==0) {
      cout << "Relative Integrated Conservation Error = ";
      for  (int i=0; i<4; i++)
       cout << std::abs(EXACT[i]-err[i])/EXACT[i] << "  ";
      cout << endl;
    }
  }

#ifndef _NO_MPI
 MPI_Finalize();
#endif
}
