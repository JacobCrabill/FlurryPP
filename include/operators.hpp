/*!
 * \file operators.hpp
 * \brief Header file for oper class
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

#include <map>
#include <vector>

#include "global.hpp"

#include "funcs.hpp"
#include "geo.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "points.hpp"

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

  /*! Setup operator to calculate divergence of correction function at solution points
   *  based upon the normal flux correction at the flux points */
  void setupCorrection(vector<point> &loc_spts, vector<point> &loc_fpts);

  //! Setup an interpolation operation between solution points and given points using solution basis
  matrix<double> setupInterpolateSptsIpts(matrix<double>& loc_ipts);

  //! Interpolate a set of values at the solution points of an element to given interpolation points
  void interpolateSptsToPoints(matrix<double>& Q_spts, matrix<double>& Q_ipts, matrix<double>& loc_ipts);

  //! Interpolate values at the solution points of an element at the given reference location
  void interpolateToPoint(matrix<double> &Q_spts, double* Q_ipts, point &loc_ipt);

  //! Interpolate a flux vector (nDims x nFields) to the given reference location
  void interpolateFluxToPoint(vector<matrix<double>> &F_spts, matrix<double> &F_ipt, point &loc_ipt);

  //! Get interpolation weights at given point
  void getBasisValues(point &ipt, vector<double> &weights);
  void getBasisValues(double* loc_ipt, double* weights);

  void setupCorrectGradU(void);

  void applyGradSpts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts);

  void applyGradFSpts(vector<matrix<double>> &F_spts, Array<matrix<double>,2>& dF_spts);

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


  void applyCorrectGradU(matrix<double>& dUc_fpts, vector<matrix<double> >& dU_spts, vector<matrix<double> > &JGinv_spts, vector<double> &detJac_spts);

  /*! Shock Capturing in the element */
  double shockCaptureInEle(matrix<double> &U_spts, double threshold);

  const matrix<double>& get_oper_div_spts();
  const matrix<double>& get_oper_spts_fpts();

  map<int,matrix<double>*> get_oper_grad_spts;
  map<int,matrix<double>*> get_oper_correct;

  //! Calculate average density over an element (needed for negative-density correction)
  void calcAvgU(matrix<double>& U_spts, vector<double>& detJ_spts, vector<double>& Uavg);

  matrix<double> interpolateCorrectedFlux(vector<matrix<double> >& F_spts, matrix<double>& dFn_fpts, point refLoc);

private:
  geo *Geo;
  input *params;
  uint nDims, nFields, eType, order, nSpts, nFpts;
  string sptsType;

  matrix<double> opp_spts_to_fpts;
  matrix<double> opp_spts_to_mpts;
  vector<matrix<double>> opp_grad_spts;
  matrix<double> opp_div_spts;
  matrix<double> opp_correction;
  vector<matrix<double>> opp_correctU;
  vector<matrix<double>> opp_correctF ;

  void setupCorrectF(vector<point> &loc_spts);

  //! Evalulate the VCJH correction function at a solution point from a flux point */
  double VCJH_quad(uint fpt, point &loc, vector<double> &spts1D, uint vcjh, uint order);
  double VCJH_hex(int fpt, point &loc, vector<double> &spts1D, uint vcjh, uint order);

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
