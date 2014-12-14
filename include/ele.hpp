/*!
 * \file ele.hpp
 * \brief Header file for ele class
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
#pragma once

#include "global.hpp"
#include "mesh.hpp"

class ele
{
public:
  int ID, IDg; //! IDg will be for MPI (if I ever get to that; for now, just a reminder!)
  int eType;
  int order;
  int nDims; // should probably just go in Geo? Or not needed once eType set?
  int nNodes;
  int nSpts, nFpts; //! # of solution points, flux points
  vector<point> loc_spts;

  mesh *Mesh; //! Pointer to mesh object to which ele 'belongs'

  ele(int in_eType, int in_order, int in_ID, vector<int> &in_nodes, mesh *in_Mesh);

  void setup(int in_eType, int in_order, int in_ID, vector<double> &xy);

  void calc_jacobian(void);

private:
  // Solution Variables
  // Still undecided on how this will be stored - double*, vector<double>, something custom?
  double* U_spts; //! Solution at solution points
  double* U_fpts; //! Solution at flux points
  double* F_spts; //! Solution at solution points
  double* F_fpts; //! Solution at flux points
  
  // Transform Variables
  double* detJac_spts;  //! Determinant of transformation Jacobian at each solution point
  double* detJac_fpts;  //! Determinant of transformation Jacobian at each solution point
  vector<vector<vector<double>>> Jac_spts;  //! Transformation Jacobian [matrix] at each solution point
  vector<vector<vector<double>>> Jac_fpts;  //! Transformation Jacobian [matrix] at each flux point
  
  // Misc.

};
