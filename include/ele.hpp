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

#include <vector>

#include "global.hpp"

//#include "face.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "matrix.hpp"

class ele
{
friend class face;
friend class boundFace;
friend class intFace;
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

  void move(int step);

  void calcGridVelocity(void);

  void calcTransforms(int initial = 0);

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

  /*! Copy U0_spts into U_spts for final time advancement */
  void copyU0_Uspts(void);
  void copyUspts_U0(void);

  /* --- Display, Output & Diagnostic Functions --- */

  /*! Get vector of primitive variables at a solution point */
  vector<double> getPrimitives(uint spt);

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

  point getPosFpt(uint spt);

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

private:

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

  // Gradients
  vector<matrix<double> > dU_spts;  //! Gradient of solution at solution points
  vector<matrix<double> > dU_fpts;  //! Gradient of solution at flux points
  vector<vector<matrix<double>>> dF_spts;  //! Gradient of flux at solution points
  vector<matrix<double>> divF_spts;         //! Divergence of flux at solution points
  vector<matrix<double>> tdF_spts;          //! Transformed gradient of flux (dF_dxi and dG_deta) at solution points

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
  vector<vector<point>> nodesRK; //! Location of mesh nodes in physical space

  // Geometry Variables
  vector<point> pos_spts;     //! Position of solution points in physical space
  vector<point> pos_fpts;     //! Position of flux points in physical space
  vector<point> pos_ppts;     //! Position of plotting points [spt+fpts+nodes]
  matrix<double> norm_fpts;   //! Unit normal in physical space
  matrix<double> tNorm_fpts;  //! Unit normal in reference space
  vector<double> dA_fpts;     //! Local equivalent face-area at flux point

  // Shock Capturing variables
  double sensor;

  // Other
  matrix<double> S_spts;      //! Entropy-adjoint variable used as error indicator for Euler
  matrix<double> S_fpts;      //! Entropy-adjoint variable at flux points
  matrix<double> S_mpts;      //! Entropy-adjoint variable at mesh points

  /* --- Temporary Variables --- */
  matrix<double> tempF;
  vector<double> tempU;

  /*! Get the values of the nodal shape bases at a solution point */
  void getShape(point loc, vector<double> &shape);

  void perturb(void);
};
