/*!
 * \file input.hpp
 * \brief Header file for input reader
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Fux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
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

  /* --- Basic Problem Variables --- */
  int equation;  //! {0 | Advection/diffusion} {1 | Euler/Navier-Stokes}
  int viscous;   //! {0 | No diffusion/Viscosity} {1 | Diffusion/Viscosity}
  int order;
  int icType;
  int motion;
  int test_case;
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
  double fix_vis;  //! Use Sutherland's Law or fixed (constant) viscosity?
  double c_sth;    //! For Sutherland's Law

  /* --- Simulation Run Parameters --- */
  int nFields;
  int nDims;
  double dt;
  double CFL;
  int dtType;
  int timeType;
  double rkTime;
  double time;
  int iterMax;
  int initIter;
  int restartIter;
  int restart;
  int restart_freq;

  int iter;

  /* --- Output Parameters --- */

  string dataFileName;
  int resType;
  int monitorResFreq;
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
  string meshFileName;
  vector<string> oversetGrids;
  int meshType;
  int nx, ny, nz;
  int nGrids;
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

  /* --- Shock Capturing Parameters --- */
  int scFlag;       // Shock Capturing Flag
  double threshold; // Threshold for considering as shock -Set to 1.0 by default

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
