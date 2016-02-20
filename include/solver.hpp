/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
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

#include <memory>
#include <map>
#include <set>
#include <vector>

class oper;
class overFace;

#include "global.hpp"

#include "ele.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "face.hpp"
#include "intFace.hpp"
#include "boundFace.hpp"
#include "mpiFace.hpp"
#include "overComm.hpp"
#include "overFace.hpp"
#include "operators.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
class tioga;
#include "tioga.h"
#endif

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
  vector<shared_ptr<ele>> eles;

  //! Vector of all non-MPI faces handled by this solver
  vector<shared_ptr<face>> faces;

  //! Vector of all MPI faces handled by this solver
  vector<shared_ptr<mpiFace>> mpiFaces;

  //! Vector of all MPI faces handled by this solver
  vector<shared_ptr<overFace>> overFaces;

  //! Local supermesh of donor elements for each cell needing to be unblanked
  vector<superMesh> donors;

#ifndef _NO_MPI
  //! Pointer to Tioga object for processing overset grids
  shared_ptr<tioga> tg;
#endif

  /* ==== Misc. Commonly-Used Variables ==== */

  int nRKSteps;

  int nGrids, gridID, gridRank, nprocPerGrid;

  int order;  //! Baseline solution order

  /* === Setup Functions === */

  solver();

  ~solver();

  //! Setup the solver with the given simulation parameters & geometry
  void setup(input *params, geo* Geo = NULL);

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
  void initializeSolution(bool PMG = false);

  void update(bool PMG_Source = false);

  //! Perform one full step of computation
  void calcResidual(int step);

  //! Calculate the stable time step limit based upon given CFL
  void calcDt(void);

  //! Advance solution in time - Generate intermediate RK stage
  void timeStepA(int step, bool PMG_Source = false);

  //! Advance solution in time - Final RK stage [assemble intermediate stages]
  void timeStepB(int step, bool PMG_Source = false);

  //! For RK time-stepping - store solution at time 'n'
  void copyUspts_U0(void);

  //! For RK time-stepping - recover solution at time 'n' before advancing to 'n+1'
  void copyU0_Uspts(void);

  //! Extrapolate the solution to the flux points
  void extrapolateU(void);

  //! Extrapolate the solution to the mesh (corner) points (and edge points in 3D)
  void extrapolateUMpts(void);

  //! Extrapolate the grid velocity to the mesh (corner) points (and edge points in 3D)
  void extrapolateGridVelMpts(void);

  //! Extrapolate entropy error-estimate variable to the mesh (corner) points (and edge points in 3D)
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

  //! Have all MPI faces begin their communication of solution data
  void doCommunication(void);

  //! Have all MPI faces begin their communication of gradient data
  void doCommunicationGrad(void);

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

  //! Calculate the viscous interface flux at all overset boundary faces
  void calcViscousFlux_overset(void);

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

  //! Integrate forces on all wall-type boundaries
  vector<double> computeWallForce(void);

  //! Integrate fluxes on all inlet/outlet boundaries [Primarily for internal flows]
  vector<double> computeMassFlux(void);

  //! For implemented test cases, calculate the L1 error over the domain
  vector<double> integrateError(void);

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

  //! For moving grids, update the overset connectivity (including adding/removing cells/faces)
  void updateOversetConnectivity(bool doBlanking = true);

  void updateOversetConnectivity2D(bool doBlanking);
  void updateOversetConnectivity3D(bool doBlanking);

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

  //! Integrate the solution error over the overset domain
  vector<double> integrateErrorOverset();

  /* ---- Callback functions specifically for TIOGA ---- */

  void getNodesPerCell(int* cellID, int* nNodes);

  void getReceptorNodes(int* cellID, int* nNodes, double* posNodes);

  void donorInclusionTest(int* cellID, double* xyz, int* passFlag, double* rst);

  void donorWeights(int* cellID, double* xyz, int* nweights, int* iNode, double* weights, double* rst, int* fracSize);

  void convertToModal(int* cellID, int* nPtsIn, double* uIn, int* nPtsOut, int* iStart, double* uOut);

  /* ---- My Overset Functions ---- */

  //! Perform Galerkin projection to fringe cells instead of boundary method
  void oversetFieldInterp();

  //! Interpolate 'right state' to overset boundaries for Riemann solve
  void oversetInterp();

  //! Interpolate 'right state' gradient to overset boundaries for viscous flux
  void oversetInterp_gradient();

  /* ---- Stabilization Functions ---- */

  void calcAvgSolution();
  bool checkDensity();
  void checkEntropy();
  void checkEntropyPlot();

  shared_ptr<overComm> OComm;

private:
  //! Pointer to the parameters object for the current solution
  input *params;

  //! Set of all element types present in current simulation
  set<int> eTypes;

  //! Set of all polynomial orders present for each element type
  map<int,set<int> > polyOrders;

  //! Lists of cells to apply various adaptation methods to
  vector<int> r_adapt_cells, h_adapt_cells, p_adapt_cells;

  /* ---- Overset Grid Variables / Functions ---- */

  vector<double> U_spts; //! Global solution vector for solver (over all elements)

  //shared_ptr<overComm> OComm;

  vector<int> iblankVert, iblankEle;
};
