/*!
 * \file flurry.hpp
 * \brief Header file for functions to run a whole simulation
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014-2016 Jacob Crabill
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#ifndef _flurry_hpp
#define _flurry_hpp

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#ifndef _NO_MPI
#include "mpi.h"
#endif

class ele;
class geo;
class solver;
class multiGrid;

#include "error.hpp"
#include "input.hpp"
#include "global.hpp"
#include "output.hpp"

class Flurry
{
public:

  /*! Assign MPI communicator and overset-grid ID
   *
   * For MPI cases, Python will call MPI_Init, so we should NOT do that here
   * Instead, Python will give us our own communicator for internal use
   * (1 communicator per grid)
   */
#ifndef _NO_MPI
  Flurry(MPI_Comm comm_in = MPI_COMM_WORLD, int n_grids = 1, int grid_id = 0);
#else
  Flurry(void);
#endif

  //! Read input file and set basic run parameters
  void read_input(char *inputfile);

  //! Perform preprocessing and prepare to run case
  void setup_solver(void);

  //! Run one full time step, including any filtering or multigrid operations
  void do_step(void);

  //! Call "do_step()" n times
  void do_n_steps(int n);

  void extrapolate_u();

  // Functions to write data to file and/or terminal
  void write_residual(void);
  void write_solution(void);
  void write_error(void);

  //void write_wall_time(void);

  // Other Misc. Functions
  input &get_input(void) { return params; }

  /* ==== Overset-Related Functions ==== */

  // Geometry Access Functions
  void get_basic_geo_data(int &btag, int &nnodes, double *&xyz, int *&iblank,
                          int &nwall, int &nover, int *&wallNodes, int *&overNodes,
                          int &nCellTypes, int &nvert_cell, int &nCells_type,
                          int *&c2v);

  void get_extra_geo_data(int &nFaceTypes, int& nvert_face, int& nFaces_type,
                          int*& f2v, int*& f2c, int*& c2f, int*& iblank_face,
                          int*& iblank_cell, int& nOver, int*& overFaces,
                          int& nMpiFaces, int*& mpiFaces, int*& procR,
                          int*& faceIdR);

  // Solution-data access functions
  double get_u_spt(int ele, int spt, int var);
  double *get_u_spts(void);
  double *get_u_fpts(void);

  double& get_u_fpt(int faceID, int fpt, int field);

  // Callback Functions for TIOGA
  void get_nodes_per_cell(int& nNodes);
  void get_nodes_per_face(int& nNodes);
  void get_receptor_nodes(int cellID, int& nNodes, double* xyz);
  void get_face_nodes(int faceID, int& nNodes, double* xyz);
  void get_q_index_face(int faceID, int fpt, int& ind, int& stride);
  void donor_inclusion_test(int cellID, double* xyz, int& passFlag, double* rst);
  void donor_frac(int cellID, int& nweights, int* inode,
                  double* weights, double* rst, int buffsize);

private:
  // Generic data about the run
  int rank = 0, nRanks = 1;
  int myGrid = 0;  //! For overset: which grid this rank belongs to
  int nGrids = 1;  //! For overset: # of grids in entire system

  // Basic ZEFR Solver Objects
  std::shared_ptr<solver> Solver;
  input params;
  std::shared_ptr<multiGrid> pmg;
  geo *Geo;

#ifndef _NO_MPI
  void mpi_init(MPI_Comm comm_in, int n_grids = 1, int grid_id = 0);
#endif

  // Again, to simplify MPI vs. no-MPI compilation, this will be either an
  // MPI_Comm or an int
  _mpi_comm myComm;
};

#endif /* _flurry_hpp */
