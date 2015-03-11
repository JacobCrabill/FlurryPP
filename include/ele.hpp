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
//#include "mesh.hpp"
#include "matrix.hpp"
#include "geo.hpp"

#include <array>

class ele
{
friend class face;
friend class bound;
friend class solver;

public:
  int ID, IDg; //! IDg will be for MPI (if I ever get to that; for now, just a reminder!)
  int eType;
  int order;
  int nNodes;

  vector<point> loc_spts; //! Location of solution points in parent domain
  vector<point> loc_fpts; //! Location of flux points in parent domain
  vector<point> nodes; //! Location of mesh nodes in physical space
  vector<int> nodeID; //! Global ID's of element's nodes
  vector<int> faceID; //! Global ID's of element's faces
  vector<bool> bndFace; //! Tag for faces on a boundary

  //! Default constructor
  ele();

  //! Alternate constructor (think I'll delete this)
  ele(int in_eType, int in_order, int in_ID, vector<point> &in_nodes, geo *in_Geo);

  void initialize(void);

  void setup(input *inParams, geo *inGeo);

  void calcTransforms(void);

  void calcPosSpts(void);

  void setInitialCondition(void);

  void calcInviscidFlux_spts(void);

  void calcViscousFlux_spts(void);

  /* --- Display, Output & Diagnostic Functions --- */

  /*! Get vector of primitive variables at a solution point */
  vector<double> getPrimitives(uint spt);

  /*! Compute the solution residual over the element */
  vector<double> getResidual(int normType);

  /*! Get position of solution point in physical space */
  point getPosSpt(uint spt);

  uint getNDims() const;
  void setNDims(int value);

  uint getNFields() const;
  void setNFields(int value);

  uint getNSpts() const;
  void setNSpts(int value);

  uint getNFpts() const;
  void setNFpts(int value);

private:

  /* --- Simulation/Mesh Parameters --- */
  geo* Geo;      //! Geometry (mesh) to which element belongs
  input* params; //! Input parameters for simulation

  int nDims;   //! # of physical dimensions for simulation
  int nFields; //! # of solution variable fields
  int nSpts;   //! # of solution points in element
  int nFpts;   //! # of flux points in element

  /* --- Solution Variables --- */
  // Solution, flux
  matrix<double> U_spts;           //! Solution at solution points
  matrix<double> U_fpts;           //! Solution at flux points
  vector<matrix<double> > F_spts;  //! Flux at solution points
  vector<matrix<double> > F_fpts;  //! Flux at flux points
  matrix<double> Fn_fpts;          //! Interface flux at flux points
  matrix<double> dFn_fpts;         //! Interface - discontinuous flux at flux points

  // Gradients
  vector<matrix<double> > dU_spts;  //! Gradient of solution at solution points
  vector<matrix<double> > dU_fpts;  //! Gradient of solution at flux points
  vector<vector<matrix<double>>> dF_spts;  //! Gradient of flux at solution points
  matrix<double> divF_spts;         //! Divergence of flux at solution points

  // Transform Variables
  vector<double> detJac_spts;  //! Determinant of transformation Jacobian at each solution point
  vector<double> detJac_fpts;  //! Determinant of transformation Jacobian at each solution point
  vector<matrix<double> > Jac_spts;  //! Transformation Jacobian [matrix] at each solution point
  vector<matrix<double> > Jac_fpts;  //! Transformation Jacobian [matrix] at each flux point
  
  // Geometry Variables
  vector<point> pos_spts;
  matrix<double> norm_fpts;   //! Unit normal in physical space
  matrix<double> tNorm_fpts;  //! Unit normal in reference space
  vector<double> dA_fpts;     //! Local equivalent face-area at flux point

  /*! Get the values of the nodal shape bases at a solution point */
  void getShape(int spt, vector<double> &shape);
};
