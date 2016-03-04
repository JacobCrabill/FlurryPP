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
  int sID;     //! ID within Solver's arrays
  int eType;
  int order;
  int nNodes; //! Number of nodes used to define element shape
  int nMpts;  //! Number of nodes used for plotting (corners, not edge nodes for quadratic eles)

  //! Default constructor
  ele();

  //! Alternate constructor (think I'll delete this)
  ele(int in_eType, int in_order, int in_ID, vector<point> &in_nodes, geo *in_Geo);

  void initialize(void);

  void setup(input *inParams, solver* inSolver, geo *inGeo, int in_order = -1);

  void move(bool doTransforms = true);

  void calcGridVelocity(void);

  void calcTransforms(bool moving = false);

  void calcTransforms_point(matrix<double>& jacobian, matrix<double>& JGinv, double& detJac, const point& loc);

  point calcPos(const point &loc);

  void setInitialCondition(void);

  void calcViscousFlux_spts(void);

  void transformGradF_spts(int step);

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
  void restart(ifstream &file, input *_params, geo *_Geo);

  void getUSpts(double* Uvec);
  void setUSpts(double* Uvec);

  bool checkDensity();
  void checkEntropy();
  void checkEntropyPlot();

  //! Calculate the error of the solution w.r.t. several test cases
  matrix<double> calcError();

  /* --- Simulation/Mesh Parameters --- */
  solver* Solver;
  geo* Geo;      //! Geometry (mesh) to which element belongs
  input* params; //! Input parameters for simulation

  int nDims;   //! # of physical dimensions for simulation
  int nFields; //! # of solution variable fields
  int nSpts;   //! # of solution points in element
  int nFpts;   //! # of flux points in element

  int nRKSteps;

  /* --- Solution Variables --- */

  vector<double> waveSp_fpts;      //! Maximum wave speed at each flux point

  vector<double> Uavg;             //! Average solution over element

  // Gradients
  vector<matrix<double>> tdF_spts;  //! Transformed gradient of flux (dF_dxi and dG_deta) at solution points

  // Shock Capturing variables
  double sensor;

  // Other
  matrix<double> S_spts;      //! Entropy-adjoint variable used as error indicator for Euler
  matrix<double> S_fpts;      //! Entropy-adjoint variable at flux points
  matrix<double> S_mpts;      //! Entropy-adjoint variable at mesh points

  /* --- Temporary Variables --- */
  double tempF[3][5];
  vector<double> tempU;

  /* --- Overset Stuff --- */
  int sptOffset;  //! Offset within overset data-transfer array to grab solution data

  vector<matrix<double>> transformFlux_physToRef(void);
  vector<matrix<double>> transformFlux_refToPhys(void);
  vector<matrix<double>> transformGradU_physToRef(void);

  /*! ==== TRANSITION TO GLOBAL SOLUTION STORAGE ==== */
  double& U_spts(int spt, int field);
  double& F_spts(int dim, int spt, int field);
  double& dU_spts(int dim, int spt, int field);
  double& dF_spts(int dim_grad, int dim_flux, int spt, int field);
  double& U_fpts(int fpt, int field);
  double& F_fpts(int dim, int fpt, int field);
  double& Fn_fpts(int fpt, int field);
  double& disFn_fpts(int fpt, int field);
  double& U_mpts(int mpt, int field);
  double& U_ppts(int ppt, int field);
  double& dU_fpts(int dim, int fpt, int field);
  double& dUc_fpts(int fpt, int field);
  double& divF_spts(int step, int spt, int field);

  double& shape_spts(int spt, int node);
  double& shape_fpts(int fpt, int node);
  double& shape_ppts(int ppt, int node);
  double& dshape_spts(int spt, int node, int dim);
  double& dshape_fpts(int fpt, int node, int dim);

  double& detJac_spts(int spt);
  double& detJac_fpts(int fpt);
  double& dA_fpts(int fpt);
  double& Jac_spts(int spt, int dim1, int dim2);
  double& Jac_fpts(int fpt, int dim1, int dim2);
  double& JGinv_spts(int spt, int dim1, int dim2);
  double& JGinv_fpts(int fpt, int dim1, int dim2);

  double& norm_fpts(int fpt, int dim);
  double& tNorm_fpts(int fpt, int dim);

  double& nodes(int npt, int dim);
  double& nodesRK(int npt, int dim);

  double& pos_spts(int spt, int dim);
  double& pos_fpts(int fpt, int dim);
  double& pos_ppts(int ppt, int dim);

  double& gridVel_spts(int spt, int dim);
  double& gridVel_fpts(int fpt, int dim);
  double& gridVel_ppts(int ppt, int dim);
  double& gridVel_nodes(int mpt, int dim);

private:

  double dt;  //! CFL-based dt for element

  /*! Get the values of the nodal shape bases at a solution point */
  void getShape(point loc, vector<double> &shape);

  double getDxNelderMead(point refLoc, point physPos);

  vector<double> tmpShape;  //! To avoid unnecessary mem allocs in calcPos
};
