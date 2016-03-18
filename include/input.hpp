/*!
 * \file input.hpp
 * \brief Header file for input reader
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 1.0.0
 *
 * Fux Reconstruction in C++ (Flurry++) Code
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "global.hpp"

class fileReader
{

public:
  /*! Default constructor */
  fileReader();

  fileReader(string fileName);

  fileReader(const fileReader& _fr);

  fileReader& operator=(const fileReader& _fr);

  /*! Default destructor */
  ~fileReader();

  /*! Set the file to be read from */
  void setFile(string fileName);

  /*! Open the file to prepare for reading simulation parameters */
  void openFile(void);

  /*! Close the file & clean up */
  void closeFile(void);

  /* === Functions to read paramters from input file === */

  /*! Read a single value from the input file; if not found, apply a default value */
  template <typename T>
  void getScalarValue(string optName, T &opt, T defaultVal);

  /*! Read a single value from the input file; if not found, throw an error and exit */
  template <typename T>
  void getScalarValue(string optName, T &opt);

  /*! Read a vector of values from the input file; if not found, apply the default value to all elements */
  template <typename T>
  void getVectorValue(string optName, vector<T> &opt, T defaultVal);

  /*! Read a vector of values from the input file; if not found, throw an error and exit */
  template <typename T>
  void getVectorValue(string optName, vector<T> &opt);

  /*! Read in a map of type <T,U> from input file; each entry prefaced by optName */
  template <typename T, typename U>
  void getMap(string optName, map<T, U> &opt);

private:
  ifstream optFile;
  string fileName;

};

class input
{
public:
  /*! Default constructor */
  input();

  void readInputFile(char *filename);

  void nonDimensionalize(void);

  simTimer timer;

  /* --- Basic Problem Variables --- */
  int equation;  //! {0 | Advection/diffusion} {1 | Euler/Navier-Stokes}
  int viscous;   //! {0 | No diffusion/Viscosity} {1 | Diffusion/Viscosity}
  int order;
  int icType;
  int motion;
  int testCase;
  int riemannType;


  /* --- Viscous Solver Parameters --- */
  double penFact;    //! Penalty factor for the LDG viscous flux
  double tau;        //! Bias parameter for the LDG viscous flux
  double Re;         //! Reynolds number
  double Lref;       //! Reference length for Reynlds number

  // For Sutherland's Law
  double muGas;
  double TGas;
  double SGas;
  double rt_inf;   //! For Sutherland's Law
  double mu_inf;   //! For Sutherland's Law
  double c_sth;    //! For Sutherland's Law
  int fixVis;  //! Use Sutherland's Law or fixed (constant) viscosity?

  /* --- Simulation Run Parameters --- */
  int nFields;
  int nDims;
  double dt;
  double CFL;
  int dtType;
  int timeType;
  double rkTime;
  double time;
  double maxTime;
  int iterMax;
  int initIter;
  int restartIter;
  int restart;
  int restart_freq;
  int nRKSteps;
  vector<double> RKa, RKb;

  /* --- Multigrid Options --- */
  int PMG;         //! P-Multigrid flag [default: off/0]
  int lowOrder;    //! Minimum order to use with PMG [default: 0]
  int smoothSteps; //! Number of 'smoothing' iterations to use on coarse levels
  int HMG;         //! H-Multigrid flag [default: off/0]
  int n_h_levels;  //! Number of h-levels to cycle [default: 1]
  int shapeOrder; //! Shape-function order to use on generated fine grids

  int iter;

  /* --- Moving-Grid Parameters --- */
  double moveAx, moveAy, moveFx, moveFy;
  double moveAr, moveFr;

  /* --- Output Parameters --- */

  string dataFileName;
  int resType;
  int monitorResFreq;
  int monitorErrFreq;
  int errorNorm;
  int quadOrder;
  int plotFreq;
  int plotType;

  bool calcEntropySensor;

  /* --- Boundary & Initial Condition Parameters --- */
  double Uinf;
  double Vinf;
  double rhoinf;
  double Pinf;

  double rhoIC;
  double vxIC, vyIC, vzIC;
  double pIC;

  double rhoBound;
  double uBound, vBound, wBound;
  double pBound;

  double TWall;

  double oneOverS;  //! Precompute for characteristic boundary condition

  // Viscous Boundary Conditions / Initial Conditions
  double nxBound, nyBound, nzBound;
  double MachBound;
  double TBound;
  double TtBound;
  double PtBound;
  double muBound;

  double TIC;
  double muIC;

  bool slipPenalty;  //! Use "penalty method" on slip-wall boundary

  /* --- Misc. Physical/Equation Parameters --- */
  double gamma;
  double prandtl;
  double mu;
  double RGas;

  double advectVx; //! Advection speed, x-direction
  double advectVy; //! Advection speed, y-direction
  double advectVz; //! Advection speed, z-direction
  double diffD;    //! Diffusion constant
  double lambda;   //! Lax-Friedrichs upwind coefficient (0: Central, 1: Upwind)

  /* --- Mesh Parameters --- */
  string meshFileName;          //! Gmsh file name for standard run
  vector<string> oversetGrids;  //! Gmsh file names of all overset grids being used
  int meshType;     //! Type of mesh being used: Single Gmsh, create a mesh, or read multiple overset grids
  int nx, ny, nz;   //! For creating a structured mesh: Number of cells in each direction
  int nGrids;       //! # of grids in overset calculation
  int writeIBLANK;  //! Write IBLANK in ParaView output?
  int oversetMethod;   //! Interp. dis. sol'n (0) or corr. flux (1) at overset bounds, or use Galerkin proj. (2) on fringe cells
  int projection;   //! Use Local Galerkin Projection (1) or simple collocation (0)
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double periodicTol, periodicDX, periodicDY, periodicDZ;
  string create_bcTop, create_bcBottom, create_bcLeft;
  string create_bcRight, create_bcFront, create_bcBack; //! BC's to apply to Flurry-created mesh
  map<string,string> meshBounds;

  /* --- FR Parameters --- */
  string sptsTypeTri;  //! Legendre, Lobatto, ...
  string sptsTypeQuad;
  int vcjhSchemeTri;
  int vcjhSchemeQuad;

  /* --- Shock Capturing, Filtering & Stabilization Parameters --- */
  int scFlag;       //! Shock Capturing Flag
  double threshold; //! Threshold for considering as shock -Set to 1.0 by default

  double exps0;     //! Minimum entropy bound for polynomial squeezing
  int squeeze;      //! Flag to turn on polynomial squeezing or not

  /* --- PID Boundary Conditions --- */
  double Kp;
  double Kd;
  double Ki;

  /* --- Other --- */
  int rank;
  int nproc;

private:
  fileReader opts;
};
