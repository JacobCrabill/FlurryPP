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

#include <sstream>

#ifndef _NO_MPI
#include <mpi.h>
#endif

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
    cout << endl;
  }

  if (argc<2) FatalError("No input file specified.");

  setGlobalVariables();

  /* Read input file & set simulation parameters */
  params.readInputFile(argv[1]);

  /* Setup the mesh and connectivity for the simulation */
  Geo.setup(&params);

  /* Setup the solver, all elements and faces, and all FR operators for computation */
  Solver.setup(&params,&Geo);

  if (!params.runOne)
    Solver.readReferenceSolution(string("cylRef"),100000);

  /* Stat timer for simulation (ignoring pre-processing) */
  simTimer runTime;
  runTime.startTimer();

  /* ---- SIMPLIFIED / MODIFIED SIMULATION RUNNING FOR AA222 FINAL PROJECT ---- */

  /* Each 'function evaluation' will entail running a standardized 2D cylinder
   * test case for 2500 iterations and recording the L2 norm of the density residual */

  params.evals = 0;
  int maxEvals = 200;

  int nVars = 3;                 // For simplex search: vars = Kp, Kd, Ki
  int nPts = nVars + 1;          // For simplex search: n_dims + 1 [# of points in nVars-dimensional simplex]
  vector<double> F(nPts);        // Objective Function Values. Using L2 norm of density after
  matrix<double> X(nPts,nVars);  // Point locations [each row is: Kp, Kd, Ki]

  params.outputPrefix = "aa222";
  params.hist.open("aa222_hist.out");
  params.hist.precision(6);
  params.hist.setf(ios::scientific, ios::floatfield);
  params.hist << "Kp, Kd, Ki, err, time" << endl;
  params.CpFile.open("Cp_hist.out");
  params.CpFile.precision(6);
  params.CpFile.setf(ios::scientific, ios::floatfield);
  params.CpFile << "Iter, Time, normPDiff" << endl;

  // Initialize the simplex from points [Kp, Kd, Ki]
  X(0,0) =  .5;  X(0,1) = 0;  X(0,2) = 0;
  X(1,0) = 1.5;  X(1,1) = 0.1;  X(1,2) = 0;
  X(2,0) =  .6;  X(2,1) = 1;  X(2,2) = 0;
  X(3,0) =  .7;  X(3,1) = .2;  X(3,2) = 1;

  cout << endl;
  cout << "AA222: Initializing Simplex with Starting Points" << endl;
  cout << endl;

  if (params.runOne) {
    params.outputPrefix = params.dataFileName;
    // For testing purposes, run once using input params and exit
    Solver.runSim(vector<double>{params.Kp,params.Kd,params.Ki},3);
    runTime.stopTimer();
    runTime.showTime();
    exit(0);
  }

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++) {
    vector<double> xTmp = {X(i,0),X(i,1),X(i,2)};
    F[i] = Solver.runSim(xTmp,3);
  }

  cout << endl;
  cout << "AA222: Running Optimization using Nelder-Meade Simplex Search" << endl;
  cout << endl;

  int iter = 0;
  while (params.evals < maxEvals) {
    cout << "Iter " << iter << endl;
    auto ind = getOrder(F);
    vector<double> Xn = {X(ind[nPts-1],0),X(ind[nPts-1],1),X(ind[nPts-1],2)};  // Point with the highest value of F
    vector<double> X0(nVars);  // Centroid of all other points
    vector<double> Xr(nVars);  // Reflected point

    for (int i=0; i<nVars; i++) {
      // Take centroid of all points besides Xn
      for (int j=0; j<nPts-1; j++) {
        X0[i] += X(ind[j],i)/(nPts-1);
      }
      // Reflect Xn around X0
      Xr[i] = X0[i] + (X0[i]-Xn[i]);
    }

    double Fr = Solver.runSim(Xr,3);

    // Determine what to do with the new point
    if (Fr < F[ind[nPts-2]]) {
      // We will be keeping this point
      if (Fr < F[ind[0]]) {
        // This one's good; keep going! Expand from Xr
        vector<double> Xe(nVars);
        for (int i=0; i<nVars; i++)
          Xe[i] = Xr[i] + (X0[i]-Xn[i]);
        double Fe = Solver.runSim(Xe,3);

        if (Fe < Fr) {
          // This one's even better; use it instead
          for (int i=0; i<nVars; i++) {
            X(ind[nPts-1],i) = Xe[i];
            F[ind[nPts-1]] = Fe;
          }
        }
        else {
          // Xe/Fe was no better; stick with Fr, Xr
          for (int i=0; i<nVars; i++) {
            X(ind[nPts-1],i) = Xr[i];
            F[ind[nPts-1]] = Fr;
          }
        }
      }
      else {
        // This one's somewhere in the middle; replace Xn with Xr
        for (int i=0; i<nVars; i++) {
          X(ind[nPts-1],i) = Xr[i];
          F[ind[nPts-1]] = Fr;
        }
      }
    }
    else {
      // Try reducing the size of the simplex
      vector<double> Xc(nVars);
      for (int i=0; i<nVars; i++) {
        Xc[i] = X0[i] - 0.5*(X0[i]-Xn[i]);
      }
      double Fc = Solver.runSim(Xc,3);
      if (Fc < F[ind[nPts-1]]) {
        // Bringing this point in is better; use it
        for (int i=0; i<nVars; i++) {
          X(ind[nPts-1],i) = Xc[i];
          F[ind[nPts-1]] = Fc;
        }
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        vector<double> X1 = {X(ind[0],0),X(ind[0],1),X(ind[0],2)};
        for (int i=1; i<nPts; i++) {
          for (int j=0; j<nVars; j++) {
            X(ind[i],j) = X1[j] + 0.5*(X(ind[i],j)-X1[j]);
          }
          vector<double> xTmp = {X(ind[i],0),X(ind[i],1),X(ind[i],2)};
          F[ind[i]] = Solver.runSim(xTmp,3);
        }
      }
    }

    // Continue to iterate
    iter++;
  }

  params.hist.close();
  params.CpFile.close();

  runTime.stopTimer();
  runTime.showTime();

#ifndef _NO_MPI
 MPI_Finalize();
#endif
}
