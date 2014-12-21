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
#include "../include/global.hpp"

#include <string>
#include <fstream>

class fileReader
{
public:
  /*! Default constructor */
  fileReader();

  /*! Default destructor */
  ~fileReader();

  /*! Open the file to prepare for reading simulation parameters */
  void openFile(string fileName);

  /*! Close the file & clean up */
  void closeFile(void);

  /* === Functions to read paramters from input file === */

  /*! Read a single value from the input file; if not found, apply a default value */
  template <typename T>
  void getScalarOpt(string optName, T &opt, T defaultVal);

  /*! Read a single value from the input file; if not found, throw an error and exit */
  template <typename T>
  void getScalarOpt(string optName, T &opt);

  /*! Read a vector of values from the input file; if not found, apply the default value to all elements */
  template <typename T>
  void getVectorVal(string optName, vector<T> &opt, T defaultVal);

  /*! Read a vector of values from the input file; if not found, throw an error and exit */
  template <typename T>
  void getVectorVal(string optName, vector<T> &opt);

private:

};

class input
{
public:
  /*! Default constructor */
  input();

  /*! Default destructor */
  ~input();

  void readInputFile(char *filename);

  /* --- Member Variables --- */
  int equation;  //! {0 | Advection/diffusion} {1 | Euler/Navier-Stokes}
  int viscous;   //! {0 | No diffusion/Viscosity} {1 | Diffusion/Viscosity}
  int order;
  int riemann_type;
  int ic_type;
  int test_case;

  /* --- Simulation Run Parameters --- */
  int nFields;
  int nDims;
  double dt;
  int iterMax;
  int plot_freq;
  int restart_freq;

  /* --- Boundary Condition Parameters --- */
  double Uinf;
  double Vinf;
  double rhoinf;
  double Pinf;

  /* --- Misc. Physical/Equation Parameters --- */
  double gamma;
  double prandtl;
  double mu;

  /* --- Mesh Parameters --- */
  string mesh_file_name;
  int mesh_type;
  int nx, ny;
  double xmin, xmax, ymin, ymax;

  /* --- FR Parameters --- */
  string spt_type_tri;  //! Legendre, Lobatto, ...
  string spt_type_quad;

private:
  fileReader opts;
};
