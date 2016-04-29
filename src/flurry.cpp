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
#include "multigrid.hpp"
#include "array"

vector<double> CycleLGPTest(Array<double,3> &U_inout, input params, string Grid1, string Grid2, int iter, double exactVal = 0);

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

#ifdef _OMP
  /* This will work, but appears to operate much less efficiently than setting
   * OMP_NUM_THREADS = # manually before run */
//  int nThreads = omp_get_num_threads();
//  omp_set_num_threads(nThreads / mpi_nproc);
#endif

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

  bool DO_LGP_TEST = true;
  if (!DO_LGP_TEST) {
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
  }

#ifndef _NO_MPI
  if (DO_LGP_TEST) {
    params.oversetGrids[0] = string("box0.msh");
    params.oversetGrids[1] = string("box1.msh");

    params.projection = 1;  //! Use LGP or direct interpolation
    params.errorNorm = 1;
    params.testCase = 0;
    params.quadOrder = 6;
    params.icType = 3;      //! Smooth Gaussian (0), sine waves (1-2), or circular step (3)

    /* Setup the solver, grid, all elements and faces, and all FR operators for computation */
    Solver.setup(&params,params.order);

    /* Apply the initial condition */
    Solver.initializeSolution();

    writeData(&Solver,&params);

    double EXACT_VAL = 0;
    if (params.icType == 0)
      EXACT_VAL = pi * erf(5) * erf(5);
    else if (params.icType == 3)
      EXACT_VAL = pi * 9;

    // Calc Error Before Interpolation:
    auto Error = Solver.integrateError();

    if (params.rank==0) EXACT_VAL = Error[0];

    MPI_Status status;
    if (params.rank == 1)
      MPI_Recv(&EXACT_VAL,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
    else
      MPI_Send(&EXACT_VAL,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);

    if (params.rank == 0) {
      cout.precision(6);
      cout.setf(ios::scientific, ios::floatfield);
      cout << endl;
      cout << "******************************************************************" << endl;
      cout << "Starting Integral Value: ";
      for (auto &val:Error) cout << val << ",  ";
      cout << endl;
      cout << "******************************************************************" << endl;
      cout << endl;
    }

    Array<double,3> U_inout;
    U_inout = Solver.U_spts;
    CycleLGPTest(U_inout, params, string("box0.msh"), string("box1.msh"), 1, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box1.msh"), string("box2.msh"), 2, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box2.msh"), string("box3.msh"), 3, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box3.msh"), string("box4.msh"), 4, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box4.msh"), string("box5.msh"), 5, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box5.msh"), string("box6.msh"), 6, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box6.msh"), string("box7.msh"), 7, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box7.msh"), string("box8.msh"), 8, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box8.msh"), string("box9.msh"), 9, EXACT_VAL);
    CycleLGPTest(U_inout, params, string("box9.msh"), string("box0.msh"), 10, EXACT_VAL);

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
  params.timer.startTimer();

  double maxTime = params.maxTime;
  int initIter = params.initIter;
  int iterMax = params.iterMax;
  int &iter = params.iter;
  iter = initIter;

  /* --- Calculation Loop --- */
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


vector<double> CycleLGPTest(Array<double,3> &U_inout, input params, string Grid1, string Grid2, int iter, double exactVal)
{
  /* --- Setup New Solver for Solution Projection from Grid1 to Grid 2 --- */

  solver Solver;
  params.oversetGrids[0] = Grid1;
  params.oversetGrids[1] = Grid2;

  Solver.setup(&params, params.order);
  geo* Geo = Solver.Geo;

  if (Geo->gridID == 0)
    Solver.U_spts = U_inout;

  Geo->fringeCells.clear();

  if (Geo->gridID==1)
    for (int ic=0; ic<Geo->nEles; ic++)
      Geo->fringeCells.insert(ic);

  // Interpolate IC to 2nd grid [from box1 to box2]
  if (params.projection)
  {
    /* --- Use Local Galerkin Projection --- */
    Solver.OComm->matchUnblankCells(Solver.eles,Geo->fringeCells,Geo->eleMap,params.quadOrder);
    Solver.OComm->performGalerkinProjection(Solver.eles,Solver.opers,Geo->eleMap,params.order);
  }
  else
  {
    /* --- Use collocation projection; Use for nonlinear shape funcs --- */
    Solver.OComm->setupFringeCellPoints(Solver.eles,Geo->fringeCells,Geo->eleMap);
    Solver.OComm->matchOversetPoints(Solver.eles,Geo->eleMap,Geo->minPt,Geo->maxPt);
    Solver.OComm->exchangeOversetData(Solver.eles,Solver.opers,Geo->eleMap);
    Solver.OComm->transferEleData(Solver.eles,Geo->fringeCells,Geo->eleMap);
  }

  params.dataFileName = "LGP_Test" + to_string(iter);
  writeData(&Solver,&params);

  // Calc Error After 1st Interpolation:
  auto Error = Solver.integrateError();

  if (params.rank == 1) {
    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    cout << endl;
    cout << "******************************************************************" << endl;
    cout << "Integrated error for iter " << iter << ": ";
    for (auto &val:Error) cout << val - exactVal << ",  ";
    cout << endl;
    cout << "******************************************************************" << endl;
    cout << endl;
  }

  vector<uint> dims = {0,0,0};
  MPI_Status status;

  if (params.rank == 0) {
    MPI_Recv(dims.data(),3,MPI_UNSIGNED,1,0,MPI_COMM_WORLD,&status);
  } else {
    U_inout = Solver.U_spts;
    dims = {U_inout.dims[0], U_inout.dims[1], U_inout.dims[2]};
    MPI_Send(dims.data(),3,MPI_UNSIGNED,0,0,MPI_COMM_WORLD);
  }

  if (params.rank == 0) {
    U_inout.setup(dims[0],dims[1],dims[2]);
    MPI_Recv(U_inout.getData(),U_inout.getSize(),MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);
  } else {
    MPI_Send(U_inout.getData(),U_inout.getSize(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }

  if (params.rank == 0) {
    MPI_Recv(Error.data(),Error.size(),MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);
  } else {
    MPI_Send(Error.data(),Error.size(),MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }

  return Error;
}
