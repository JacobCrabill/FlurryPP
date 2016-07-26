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
#include <chrono>

#ifndef _NO_MPI
#include <mpi.h>
#endif

#ifndef _NO_MPI_DEBUG
#include <unistd.h>  // for getpid()
#endif

#include "funcs.hpp"
#include "input.hpp"
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

#ifndef _NO_MPI_DEBUG
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

  params.time1.setPrefix("MPI Wait Time: ");
  params.time2.setPrefix("Overset MPI Time: ");

  CALLGRIND_START_INSTRUMENTATION;
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
  CALLGRIND_STOP_INSTRUMENTATION;
  /* Calculate the integral / L1 / L2 error for the final time */
  writeAllError(&Solver,&params);
MPI_Barrier(MPI_COMM_WORLD);

//if (rank == 0)
//{
  std::cout << "Interp Time: ";
  params.interpTime.showTime(2);
  std::cout << "Computation Time: ";
  params.runTime.showTime(2);
//}

  // Get simulation wall time
  params.timer.stopTimer();
  params.timer.showTime();

  params.time1.showTime(6);
  params.time2.showTime(6);

#ifndef _NO_MPI
 MPI_Finalize();
#endif

 return 0;
}



/* ==== Add in interface functions for use from external code ==== */
//#define _BUILD_LIB  // for QT editing
//#ifdef _BUILD_LIB

#ifndef _NO_MPI
Flurry::Flurry(MPI_Comm comm_in, int n_grids, int grid_id)
#else
Flurry::Flurry(void)
#endif
{
  // Basic constructor
#ifdef _NO_MPI
  myComm = 0;
  myGrid = 0;
#else
  mpi_init(comm_in, n_grids, grid_id);
#endif

  /* Print out cool ascii art header */
  if (myGrid == 0 && rank == 0)
  {
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

  Solver = NULL;
}

#ifndef _NO_MPI
void Flurry::mpi_init(MPI_Comm comm_in, int n_grids, int grid_id)
{
  myComm = comm_in;
  nGrids = n_grids;
  myGrid = grid_id;

  MPI_Comm_rank(myComm, &rank);
  MPI_Comm_size(myComm, &nRanks);

#ifdef _GPU
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  /// TODO
  if (nDevices < nRanks)
  {
    //ThrowException("Not enough GPUs for this run. Allocate more!");
  }

  cudaSetDevice(rank%6); // Hardcoded for ICME nodes for now.
#endif
}
#endif

void Flurry::read_input(char *inputfile)
{
  if (rank == 0) std::cout << "Reading input file: " << inputfile <<  std::endl;
  params.readInputFile(inputfile);

  params.rank = rank;
  params.nproc = nRanks;
  params.gridID = myGrid;
  params.nGrids = nGrids;
  params.meshType = READ_MESH;
  params.myComm = myComm;
  if (nGrids > 1) params.overset = true;

  if (params.overset)
  {
    params.meshFileName = params.oversetGrids[myGrid];
  }
}

void Flurry::setup_solver(void)
{
  if (!Solver)
    Solver = make_shared<solver>();

  Solver->myComm = myComm;
  Solver->gridID = myGrid;

  if (params.PMG)
  {
    /* Setup the P-Multigrid class if requested */
    pmg = make_shared<multiGrid>();
    pmg->setup(params.order,&params,*Solver);
  }
  else
  {
    /* Setup the solver, grid, all elements and faces, and all FR operators for computation */
    Solver->setup(&params,params.order);

    /* Apply the initial condition */
    Solver->initializeSolution();
  }

  Geo = Solver->Geo;
}

void Flurry::do_step(void)
{
  params.iter++;

  params.timer.startTimer();

  Solver->update();

  /* If using multigrid, perform correction cycle */
  if (params.PMG)
    pmg->cycle(*Solver);

  params.timer.stopTimer();
}

void Flurry::do_n_steps(int n)
{
  for (int i = 0; i < n; i++)
    do_step();
}

void Flurry::extrapolate_u(void)
{
  Solver->extrapolateU();
}

void Flurry::write_residual(void)
{
  writeResidual(Solver.get(), &params);
}

void Flurry::write_solution(void)
{
  writeData(Solver.get(), &params);
}

void Flurry::write_error(void)
{
  writeError(Solver.get(), &params);
}

void Flurry::get_basic_geo_data(int& btag, int& nnodes, double*& xyz, int*& iblank,
                              int& nwall, int& nover, int*& wallNodes,
                              int*& overNodes, int& nCellTypes, int &nvert_cell,
                              int &nCells_type, int*& c2v)
{
  btag = myGrid;
  nnodes = Geo->nVerts;
  xyz = Geo->xv.getData();
  iblank = Geo->iblank.data();
  nwall = Geo->iwall.size();
  nover = Geo->iover.size();
  wallNodes = Geo->iwall.data();
  overNodes = Geo->iover.data();
  nCellTypes = 1;
  nvert_cell = *Geo->nodesPerCell;
  nCells_type = Geo->nEles;
  c2v = Geo->c2v.getData();
}

void Flurry::get_extra_geo_data(int& nFaceTypes, int& nvert_face,
                              int& nFaces_type, int*& f2v, int*& f2c, int*& c2f,
                              int*& iblank_face, int*& iblank_cell,
                              int &nOver, int*& overFaces, int &nMpiFaces, int*& mpiFaces, int*& procR,
                              int*& faceIdR)
{
  nFaceTypes = 1;
  nvert_face = Geo->f2nv[0];
  nFaces_type = Geo->nFaces;
  f2v = Geo->f2v.getData();
  f2c = Geo->f2c.getData();
  c2f = Geo->c2f.getData();
  iblank_face = Geo->iblankFace.data();
  iblank_cell = Geo->iblankCell.data();
  nOver = Geo->fringeFaces.size();
  overFaces = Geo->fringeFaces.data();
  nMpiFaces = Geo->nMpiFaces;
  mpiFaces = Geo->mpiFaces.data();
  procR = Geo->procR.data();
  faceIdR = Geo->faceID_R.data();
}

double Flurry::get_u_spt(int ele, int spt, int var)
{
  return Solver->U_spts(spt, ele, var);
}

double *Flurry::get_u_spts(void)
{
  return Solver->U_spts.getData();
}

double *Flurry::get_u_fpts(void)
{
  return Solver->U_fpts.getData();
}

void Flurry::get_nodes_per_cell(int &nNodes)
{
  nNodes = (int)Solver->nSpts;
}

void Flurry::get_nodes_per_face(int& nNodes)
{
  nNodes = (int) (Solver->nFpts / Geo->c2nf[0]);
}

void Flurry::get_receptor_nodes(int cellID, int& nNodes, double* xyz)
{
  nNodes = (int)Solver->nSpts;

  for (int spt = 0; spt < nNodes; spt++)
    for (int dim = 0; dim < Geo->nDims; dim++)
      xyz[3*spt+dim] = Solver->pos_spts(spt, cellID, dim);
}

void Flurry::get_face_nodes(int faceID, int &nNodes, double* xyz)
{
  nNodes = (int)(Solver->nFpts / Geo->c2nf[0]);

  int ff = Geo->faceMap[faceID];

  for (int fpt = 0; fpt < nNodes; fpt++)
  {
    point pt;
    if (Geo->faceType[faceID] == MPI_FACE)
      pt = Solver->mpiFaces[ff]->getPosFpt(fpt);
    else
      pt = Solver->faces[ff]->getPosFpt(fpt);

    for (int dim = 0; dim < Geo->nDims; dim++)
      xyz[3*fpt+dim] = pt[dim];
  }
}

void Flurry::get_q_index_face(int faceID, int fpt, int& ind, int& stride)
{
  int ff = Geo->faceMap[faceID];
//cout << "rank " << params.rank << " ff = " << ff << endl;

  if (Geo->faceType[faceID] == MPI_FACE)
  {
//    cout << "nMpiFaces = " << Solver->mpiFaces.size() << endl;
    Solver->mpiFaces[ff]->get_U_index(fpt,ind,stride);
  }else{
//    cout << "nFaces = " << Solver->faces.size() << endl;
    Solver->faces[ff]->get_U_index(fpt,ind,stride);
  }
}

double& Flurry::get_u_fpt(int faceID, int fpt, int field)
{
  int ff = Geo->faceMap[faceID];
  if (Geo->faceType[faceID] == MPI_FACE)
    return Solver->mpiFaces[ff]->get_U_fpt(fpt,field);
  else
    return Solver->faces[ff]->get_U_fpt(fpt,field);
}

void Flurry::donor_inclusion_test(int cellID, double* xyz, int& passFlag, double* rst)
{
  point pos = point(xyz);
  point loc;

  passFlag = Solver->eles[cellID]->getRefLocNewton(pos,loc);

  rst[0] = loc.x;
  rst[1] = loc.y;
  rst[2] = loc.z;
}

void Flurry::donor_frac(int cellID, int &nweights, int* inode, double* weights,
                      double* rst, int buffsize)
{
  Solver->opers[Solver->eles[cellID]->order].getBasisValues(rst, weights);
  nweights = Solver->eles[cellID]->nSpts;
}

//#endif /* _BUILD_LIB */
