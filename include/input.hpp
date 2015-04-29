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

  /* --- Basic Problem Variables --- */
  int equation;  //! {0 | Advection/diffusion} {1 | Euler/Navier-Stokes}
  int viscous;   //! {0 | No diffusion/Viscosity} {1 | Diffusion/Viscosity}
  int order;
  int riemann_type;
  int ic_type;
  int test_case;
  int motion;

  /* --- Simulation Run Parameters --- */
  int nFields;
  int nDims;
  double dt;
  int iterMax;
  int initIter;
  int restartIter;
  int restart;
  int plot_freq;
  int restart_freq;
  int resType;
  int monitor_res_freq;

  string dataFileName;

  /* --- Boundary & Initial Condition Parameters --- */
  double Uinf;
  double Vinf;
  double rhoinf;
  double Pinf;

  double rhoIC;
  double vxIC;
  double vyIC;
  double pIC;

  double rhoBound, uBound, vBound, wBound, pBound, TBound;
  double TWall;

  /* --- Misc. Physical/Equation Parameters --- */
  double gamma;
  double prandtl;
  double mu;
  double RGas;

  double rt_inf;  // ?
  double mu_inf;
  double fix_vis;
  double c_sth;

  double advectVx;
  double advectVy;
  double lambda;  //! Lax-Friedrichs upwind coefficient (0: Central, 1: Upwind)

  /* --- Mesh Parameters --- */
  string meshFileName;
  int mesh_type;
  int nx, ny;
  double xmin, xmax, ymin, ymax;
  double periodicTol, periodicDX, periodicDY;
  string create_bcTop, create_bcBottom, create_bcLeft, create_bcRight; //! BC's to apply to Flurry-created mesh
  //map<string,int> bcNum;
  map<string,string> meshBounds;

  /* --- FR Parameters --- */
  string sptsTypeTri;  //! Legendre, Lobatto, ...
  string sptsTypeQuad;
  int vcjhSchemeTri;
  int vcjhSchemeQuad;

private:
  fileReader opts;
};
