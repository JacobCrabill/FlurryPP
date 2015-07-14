/*!
 * \file operators.hpp
 * \brief Header file for oper class
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

#include <map>
#include <vector>

#include "global.hpp"

#include "geo.hpp"
#include "input.hpp"
#include "matrix.hpp"

class oper
{
public:
  //! Overall setup function for one element type & polynomial order
  void setupOperators(uint eType, uint order, geo* inGeo, input* inParams);

  //! Setup operator for extrapolation from solution points to flux points
  void setupExtrapolateSptsFpts(vector<point> &loc_fpts);

  //! Setup operator for extrapolation from solution points to mesh (corner) points
  void setupExtrapolateSptsMpts(vector<point> &loc_spts);

  //! Setup operator for calculation of gradient at the solution points
  void setupGradSpts(vector<point> &loc_spts);

  //! Setup an interpolation operation between two sets of points using solution basis
  void setupInterpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp);

  /*! Setup operator to calculate divergence of correction function at solution points
   *  based upon the normal flux correction at the flux points */
  void setupCorrection(vector<point> &loc_spts, vector<point> &loc_fpts);


  void setupCorrectGradU(void);

  void applyGradSpts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts);

  void applyGradFSpts(vector<matrix<double>> &F_spts, vector<vector<matrix<double>>> &dF_spts);

  void applyDivFSpts(vector<matrix<double>> &F_spts, matrix<double> &divF_spts);

  void applySptsFpts(matrix<double> &U_spts, matrix<double> &U_fpts);

  void applySptsMpts(matrix<double> &U_spts, matrix<double> &U_mpts);

  /*! For the standard FR method: extrapolate the transformed flux to the flux points
   *  and dot with the transformed outward unit normal */
  void applyExtrapolateFn(vector<matrix<double>> &F_spts, matrix<double> &tnorm_fpts, matrix<double> &Fn_fpts);

  /*! For the modified space-time transformation method: Extrapolate the physical flux
   *  to the flux points and dott with the physical outward unit normal */
  void applyExtrapolateFn(vector<matrix<double>> &F_spts, matrix<double> &norm_fpts, matrix<double> &Fn_fpts, vector<double> &dA_fpts);


  void applyCorrectDivF(matrix<double> &dFn_fpts, matrix<double> &divF_spts);


  void applyCorrectGradU(matrix<double>& dUc_fpts, vector<matrix<double> >& dU_spts);

  /*! Shock Capturing in the element */
  double shockCaptureInEle(matrix<double> &U_spts, double threshold);

  const matrix<double>& get_oper_div_spts();
  const matrix<double>& get_oper_spts_fpts();

  map<int,matrix<double>*> get_oper_grad_spts;
  map<int,matrix<double>*> get_oper_correct;

  //! Calculate average density over an element (needed for negative-density correction)
  void calcAvgU(matrix<double>& U_spts, vector<double>& detJ_spts, vector<double>& Uavg);

private:
  geo *Geo;
  input *params;
  uint nDims, nFields, eType, order, nSpts, nFpts;

  matrix<double> opp_spts_to_fpts;
  matrix<double> opp_spts_to_mpts;
  vector<matrix<double>> opp_grad_spts;
  matrix<double> opp_div_spts;
  matrix<double> opp_correction;
  vector<matrix<double>> opp_correctU;

  /*! Evaluate the divergence of the (VCJH) correction function at a solution point from a flux point */
  double divVCJH_quad(int in_fpt, point& loc, vector<double> &loc_1d_spts, uint vcjh, uint order);

  double divVCJH_hex(int in_fpt, point& loc, vector<double> &loc_1d_spts, uint vcjh, uint order);

  /* Stuff required for Shock Capturing */
  matrix<double> vandermonde1D;
  matrix<double> inv_vandermonde1D;
  matrix<double> vandermonde2D;
  matrix<double> inv_vandermonde2D;
  matrix<double> sensingMatrix;
  matrix<double> filterMatrix;
  void setupVandermonde(vector<point> &loc_spts);
  void setupSensingMatrix(void);
  void setupFilterMatrix(void);
};
