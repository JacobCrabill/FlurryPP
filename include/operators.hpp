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

#include "global.hpp"
//#include "mesh.hpp"
#include "geo.hpp"
#include "polynomials.hpp"

class oper
{
public:
  //! Overall setup function for one element type & polynomial order
  void setupOperators(uint eType, uint order, geo* inGeo, input* inParams);

  //! Setup operator for extrapolation from solution points to flux points
  void setupExtrapolateSptsFpts(vector<point> loc_spts, vector<point> loc_fpts);

  //! Setup operator for calculation of gradient at the solution points
  void setupGradSpts(vector<point> loc_spts);

  //! Setup an interpolation operation between two sets of points using solution basis
  void setupInterpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp);

  /*! Setup operator to calculate divergence of correction function at solution points
   *  based upon the normal flux correction at the flux points */
  void setupCorrection(vector<point> loc_spts, vector<point> loc_fpts);

  // Create a map<int,double*> (?) to get access to the correct operator
  // i.e. somthing like: div_flux_spts_tri = oper.get_oper_div[TRI]

  void applyGradSpts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts);

  void applyGradFSpts(vector<matrix<double>> &F_spts, vector<vector<matrix<double>>> &dF_spts);

  void applyDivFSpts(vector<matrix<double>> &F_spts, matrix<double> &divF_spts);

  void applySptsFpts(matrix<double> &U_spts, matrix<double> &U_fpts);

  void applyExtrapolateFn(vector<matrix<double>> &F_spts, matrix<double> &norm_fpts, matrix<double> &Fn_fpts);

  void applyCorrectDivF(matrix<double> &dFn_fpts, matrix<double> &divF_spts);

  const matrix<double>& get_oper_div_spts();
  const matrix<double>& get_oper_spts_fpts();

  //map<int,matrix<double>*> get_oper_div_spts;
  map<int,matrix<double>*> get_oper_grad_spts;
  //map<int,matrix<double>*> get_oper_spts_fpts;
  map<int,matrix<double>*> get_oper_correct;

private:
  geo *Geo;
  input *params;
  uint nDims, nFields, eType, order;

  matrix<double> opp_spts_to_fpts;
  vector<matrix<double>> opp_grad_spts;
  matrix<double> opp_div_spts;
  matrix<double> opp_correction;

  /*! Evaluate the divergence of the (VCJH) correction function at a solution point from a flux point */
  double divVCJH_quad(int in_fpt, vector<double> &loc, vector<double> &loc_1d_spts, uint vcjh, uint order);
};
