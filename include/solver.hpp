/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include <memory>
#include <map>
#include <set>
#include <vector>

class oper;

#include "global.hpp"

#include "ele.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "face.hpp"
#include "intFace.hpp"
#include "boundFace.hpp"
#include "mpiFace.hpp"
#include "overFace.hpp"
#include "operators.hpp"

class tioga;

#include "tioga.h"

class solver
{
friend class geo; // Probably only needed if I make eles, opers private?
friend class overFace;

public:
  /* === Member Variables === */
  //! Pointer to geometry object for mesh-related operations
  geo *Geo;

  //! Map from eType to order to element- & order-specific operator
  map<int, map<int,oper> > opers;

  //! Vector of all eles handled by this solver
  vector<ele> eles;

  //! Vector of all non-MPI faces handled by this solver
  vector<shared_ptr<face>> faces;

  //! Vector of all MPI faces handled by this solver
  vector<shared_ptr<mpiFace>> mpiFaces;

  //! Vector of all MPI faces handled by this solver
  vector<shared_ptr<overFace>> overFaces;

  //! Pointer to Tioga object for processing overset grids
  tioga* tg;

  /* === Setup Functions === */

  solver();

  ~solver();

  //! Setup the solver with the given simulation parameters & geometry
  void setup(input *params, geo *Geo);

  //! Setup the FR operators for all ele types and polynomial orders which will be used in computation
  void setupOperators();

  //! Run the basic setup functions for all elements and faces
  void setupElesFaces();

  //! If restarting from data file, read data and setup eles & faces accordingly
  void readRestartFile();

  //! Finish setting up the MPI faces
  void finishMpiSetup(void);

  /* === Functions Related to Basic FR Process === */

  //! Apply the initial condition to all elements
  void initializeSolution();

  void update(void);

  //! Perform one full step of computation
  void calcResidual(int step);

  //! Calculate the stable time step limit based upon given CFL
  void calcDt(void);

  //! Advance solution in time - Generate intermediate RK stage
  void timeStepA(int step);

  //! Advance solution in time - Final RK stage [assemble intermediate stages]
  void timeStepB(int step);

  //! For RK time-stepping - store solution at time 'n'
  void copyUspts_U0(void);

  //! For RK time-stepping - recover solution at time 'n' before advancing to 'n+1'
  void copyU0_Uspts(void);

  //! Extrapolate the solution to the flux points
  void extrapolateU(void);

  //! Extrapolate the solution to the mesh (corner) points
  void extrapolateUMpts(void);

  //! Extrapolate entropy error-estimate variable to the mesh (corner) points
  void extrapolateSMpts(void);

  //! //! Extrapolate entropy error-estimate variable to the flux points
  void extrapolateSFpts(void);

  //! Calculate the inviscid flux at the solution points
  void calcInviscidFlux_spts(void);

  //! Calculate the inviscid interface flux at all element faces
  void calcInviscidFlux_faces(void);

  //! Calculate the inviscid interface flux at all MPI-boundary faces
  void calcInviscidFlux_mpi(void);

  //! Calculate the inviscid interface flux at all overset-boundary faces
  void calcInviscidFlux_overset(void);

  //! Have all MPI faces begin their communication
  void doCommunication(void);

  //! Calculate the gradient of the solution at the solution points
  void calcGradU_spts(void);

  //! For viscous calculations, apply the correction procedure to the solution
  void correctGradU(void);

  /*! For viscous calculations, extrapolate the corrected gradient of the solution
   *  from the solution points to the flux points */
  void extrapolateGradU(void);

  //! Calculate the viscous flux at the solution points
  void calcViscousFlux_spts(void);

  //! Calculate the viscous interface flux at all element faces
  void calcViscousFlux_faces(void);

  //! Calculate the viscous interface flux at all MPI boundary faces
  void calcViscousFlux_mpi(void);

  //! Wrapper to calc divergence of flux (using one of two possible methods)
  void calcFluxDivergence(int step);

  //! Calculate the gradient of the flux at the solution points
  void calcGradF_spts(void);

  //! Use full space-time chain-rule to transform gradients to parent domain
  void transformGradF_spts(int step);

  //! Calculate the divergence of the flux at the solution points
  void calcDivF_spts(int step);

  //! Extrapolate total flux to flux points & dot with normal
  void extrapolateNormalFlux(void);

  //! Apply the correction function & add to the divergence of the flux
  void correctDivFlux(int step);

  //! Apply mesh motion
  void moveMesh(int step);

  vector<double> computeWallForce(void);

  /* === Functions for Shock Capturing & Filtering=== */

  //! Use concentration sensor + exponential modal filter to capture discontinuities
  void shockCapture(void);

  /* === Functions Related to Adaptation === */

  //! Calculate an entropy-adjoint-based error indicator
  void calcEntropyErr_spts();

  // **All of the following functions are just food for thought at the moment**

  void get_r_adapt_cells();

  void get_p_adapt_cells();

  void get_h_adapt_cells();

  void setup_r_adaption();

  void setup_h_adaptation();

  void setup_p_adaptation();

  void add_ele(int eType, int order);

  /* === Functions Related to Overset Grids === */

  //! Do initial preprocessing for overset grids
  void setupOverset();

  /*!
   * \brief Initialize overset-related data storage
   *
   * Allocates global storage for overset data. Re-call if adding elements,
   * changing polynomial orders, etc.
   */
  void setupOversetData();

  //! Copy data from elements into global solution array
  void setGlobalSolutionArray();

  //! Copy data from global solution array back into elements
  void updateElesSolutionArrays();

  //! Perform the overset data interpolation using TIOGA high-order
  void callDataUpdateTIOGA();

  /* ---- Callback functions specifically for TIOGA ---- */

  void getNodesPerCell(int* cellID, int* nNodes);

  void getReceptorNodes(int* cellID, int* nNodes, double* posNodes);

  void donorInclusionTest(int* cellID, double* xyz, int* passFlag, double* rst);

  void donorWeights(int* cellID, double* xyz, int* nweights, int* iNode, double* weights, double* rst, int* fracSize);

  void convertToModal(int* cellID, int* nPtsIn, double* uIn, int* nPtsOut, int* iStart, double* uOut);

  /* ---- My Overset Functions ---- */

  void oversetInterp();
  void sendRecvOversetData();

  /* ---- Stabilization Functions ---- */

  void calcAvgSolution();
  bool checkDensity();
  void checkEntropy();
  void checkEntropyPlot();

private:
  //! Pointer to the parameters object for the current solution
  input *params;

  //! Set of all element types present in current simulation
  set<int> eTypes;

  //! Set of all polynomial orders present for each element type
  map<int,set<int> > polyOrders;

  //! Lists of cells to apply various adaptation methods to
  vector<int> r_adapt_cells, h_adapt_cells, p_adapt_cells;

  int nRKSteps;

  int gridID, gridRank, nprocPerGrid;

  vector<double> RKa, RKb;

  /* ---- Overset Grid Variables / Functions ---- */

  vector<double> U_spts; //! Global solution vector for solver (over all elements)

  // Outgoing (interpolated) data
  vector<int> interpProc;       //! Processor ID for all interpolation points
  vector<int> interpCell;       //! For each interpolation point, cell which it lies within
  matrix<double> interpPtsPhys; //! Physical positions of fringe points on other grids to interpolate to
  matrix<double> interpPtsRef;  //! Reference position (within interpCell) of fringe points on other grids to interpolate to
  vector<matrix<double>> U_ipts; //! Interpolated solution vector at each interpolation point for each grid

  // Incoming (overset) data
  vector<int> overProc;         //! Donor-processor ID for each fringe point
  matrix<double> overPtsPhys;   //! Physical positions of each fringe point
  matrix<double> U_opts;        //! Solution vector at each fringe point
  vector<point> overPts;        //! Physical positions of fringe points
  int nOverPts;                 //! Number of overset flux points on this grid rank
};
