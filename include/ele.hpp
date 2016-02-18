/*!
 * \file ele.hpp
 * \brief Header file for ele class
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include <vector>

#include "global.hpp"

#include "funcs.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "points.hpp"

class ele
{
friend class face;
friend class boundFace;
friend class intFace;
friend class overFace;
friend class solver;

public:
  int ID, IDg; //! IDg is global ID in MPI cases
  int eType;
  int order;
  int nNodes; //! Number of nodes used to define element shape
  int nMpts;  //! Number of nodes used for plotting (corners, not edge nodes for quadratic eles)

  string sptsType;  //! Which set of point locations to use for solution and flux points

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

  void move(bool doTransforms = true);

  void calcGridVelocity(void);

  void calcTransforms(bool moving = false);

  void calcTransforms_point(matrix<double>& jacobian, matrix<double>& JGinv, double& detJac, const point& loc);

  point calcPos(const point &loc);

  void calcPosSpts(void);

  void calcPosFpts(void);

  void updatePosSpts(void);

  void updatePosFpts(void);

  void setPpts(void);

  void setShape_spts(void);

  void setShape_fpts(void);

  void setDShape_spts(void);

  void setDShape_fpts(void);

  void setTransformedNormals_fpts(void);

  void setInitialCondition(void);

  void calcInviscidFlux_spts(void);

  void calcViscousFlux_spts(void);

  void transformGradF_spts(int step);

  void calcDeltaFn(void);

  void calcDeltaUc(void);

  /*! Calculate the maximum stable time step based upon CFL */
  double calcDt(void);

  /*! Calculate the wave speed for use with calculating allowable DT */
  void calcWaveSpFpts();

  /*! Advance intermediate stages of Runge-Kutta time integration */
  void timeStepA(int step, double rkVal);

  /*! Perform final advancement of Runge-Kutta time integration */
  void timeStepB(int step, double rkVal);

  /*! Advance intermediate stages of Runge-Kutta time integration [With PMG source term] */
  void timeStepA_source(int step, double rkVal);

  /*! Perform final advancement of Runge-Kutta time integration [With PMG source term] */
  void timeStepB_source(int step, double rkVal);

  /*! Copy U0_spts into U_spts for final time advancement */
  void copyU0_Uspts(void);
  void copyUspts_U0(void);

  /* --- Display, Output & Diagnostic Functions --- */

  /*! Get vector of primitive variables at a solution point */
  vector<double> getPrimitives(uint spt);

  /*! Get vector of primitive variables at a flux point */
  vector<double> getPrimitivesFpt(uint fpt);

  /*! Get vector of primitive variables at a mesh point */
  vector<double> getPrimitivesMpt(uint mpt);

  /*! Get the full matrix of solution values at spts + fpts combined */
  void getPrimitivesPlot(matrix<double> &V);

  /*! Get the full set of grid velocity values at spts + fpts combined */
  void getGridVelPlot(matrix<double> &GV);

  /*! Get the locations of the plotting points */
  vector<point> getPpts(void);

  /*! Compute the norm of the solution residual over the element */
  vector<double> getNormResidual(int normType);

  /*! Get position of solution point in physical space */
  point getPosSpt(uint spt);

  point getPosFpt(uint fpt);

  void getPosSpts(double* posSpts);

  vector<point> getPosSpts();

  //! Get a bounding box for the element defined by the minimum and maximum extents
  vector<double> getBoundingBox(void);

  /*! Find the reference location of a point inside an element given its
   *  physical location, using the Newton root-finding method */
  bool getRefLocNewton(point pos, point& loc);

  /*! Find the reference location of a point inside an element given its
   *  physical location, using the Nelder-Meade algorithm */
  bool getRefLocNelderMead(point pos, point &loc);

  uint getNDims() const;
  void setNDims(int value);

  uint getNFields() const;
  void setNFields(int value);

  uint getNSpts() const;
  void setNSpts(int value);

  uint getNFpts() const;
  void setNFpts(int value);

  double getSensor(void);

  void calcEntropyErr_spts(void);
  vector<double> getEntropyVars(int spt);
  void getEntropyErrPlot(matrix<double> &S);
  void setupArrays();
  void setupAllGeometry();
  void restart(ifstream &file, input *_params, geo *_Geo);

  void getUSpts(double* Uvec);
  void setUSpts(double* Uvec);

  bool checkDensity();
  void checkEntropy();
  void checkEntropyPlot();

  //! Calculate the error of the solution w.r.t. several test cases
  matrix<double> calcError();

  /* --- Simulation/Mesh Parameters --- */
  geo* Geo;      //! Geometry (mesh) to which element belongs
  input* params; //! Input parameters for simulation

  int nDims;   //! # of physical dimensions for simulation
  int nFields; //! # of solution variable fields
  int nSpts;   //! # of solution points in element
  int nFpts;   //! # of flux points in element

  int nRKSteps;

  /* --- Solution Variables --- */
  // Solution, flux
  matrix<double> U_spts;           //! Solution at solution points
  matrix<double> U_fpts;           //! Solution at flux points
  matrix<double> U_mpts;           //! Solution at mesh (corner) points
  matrix<double> U0;               //! Solution at solution points, beginning of each time step
  vector<matrix<double> > F_spts;  //! Flux at solution points
  vector<matrix<double> > F_fpts;  //! Flux at flux points
  matrix<double> disFn_fpts;       //! Discontinuous normal flux at flux points
  matrix<double> Fn_fpts;          //! Interface flux at flux points
  matrix<double> dFn_fpts;         //! Interface minus discontinuous flux at flux points
  matrix<double> Uc_fpts;          //! Common solution at flux points
  matrix<double> dUc_fpts;         //! Common minus discontinuous solution at flux points
  vector<double> waveSp_fpts;      //! Maximum wave speed at each flux point

  vector<double> Uavg;             //! Average solution over element

  // Gradients
  vector<matrix<double> > dU_spts;  //! Gradient of solution at solution points
  vector<matrix<double> > dU_fpts;  //! Gradient of solution at flux points
  Array<matrix<double>,2> dF_spts;  //! Gradient of flux at solution points
  vector<matrix<double>> divF_spts; //! Divergence of flux at solution points
  vector<matrix<double>> tdF_spts;  //! Transformed gradient of flux (dF_dxi and dG_deta) at solution points

  // Transform Variables
  vector<double> detJac_spts;  //! Determinant of transformation Jacobian at each solution point
  vector<double> detJac_fpts;  //! Determinant of transformation Jacobian at each solution point
  vector<matrix<double> > Jac_spts;  //! Transformation Jacobian [matrix] at each solution point
  vector<matrix<double> > Jac_fpts;  //! Transformation Jacobian [matrix] at each flux point
  vector<matrix<double> > JGinv_spts;  //! Inverse of transformation Jacobian [matrix] at each solution point
  vector<matrix<double> > JGinv_fpts;  //! Inverse of transformation Jacobian [matrix] at each flux point

  matrix<double> shape_spts;
  matrix<double> shape_fpts;
  vector<matrix<double>> dShape_spts;  //! Derivative of shape basis at solution points
  vector<matrix<double>> dShape_fpts;  //! Derivative of shape basis at flux points
  matrix<double> gridVel_spts;         //! Mesh velocity at solution points
  matrix<double> gridVel_fpts;         //! Mesh velocity at flux points
  matrix<double> gridVel_nodes;        //! Mesh velocity at mesh (corner) points
  matrix<double> gridVel_mpts;         //! Mesh velocity at ALL mesh points (corners + edges in 3D)
  vector<point> nodesRK;               //! Location of mesh nodes in physical space

  // Geometry Variables
  vector<point> pos_spts;     //! Position of solution points in physical space
  vector<point> pos_fpts;     //! Position of flux points in physical space
  vector<point> pos_ppts;     //! Position of plotting points [spt+fpts+nodes]
  matrix<double> norm_fpts;   //! Unit normal in physical space
  matrix<double> tNorm_fpts;  //! Unit normal in reference space
  vector<double> dA_fpts;     //! Local equivalent face-area at flux point

  // Shock Capturing variables
  double sensor;

  // Multigrid Variables
  matrix<double> corr_spts, src_spts, sol_spts;

  // Other
  matrix<double> S_spts;      //! Entropy-adjoint variable used as error indicator for Euler
  matrix<double> S_fpts;      //! Entropy-adjoint variable at flux points
  matrix<double> S_mpts;      //! Entropy-adjoint variable at mesh points

  /* --- Temporary Variables --- */
  matrix<double> tempF;
  vector<double> tempU;

  /* --- Overset Stuff --- */
  int sptOffset;  //! Offset within overset data-transfer array to grab solution data

  vector<matrix<double>> transformFlux_physToRef(void);
  vector<matrix<double>> transformFlux_refToPhys(void);
  vector<matrix<double>> transformGradU_physToRef(void);

private:

  double dt;  //! CFL-based dt for element

  /*! Get the values of the nodal shape bases at a solution point */
  void getShape(point loc, vector<double> &shape);

  double getDxNelderMead(point refLoc, point physPos);

  vector<double> tmpShape;  //! To avoid unnecessary mem allocs in calcPos
};
