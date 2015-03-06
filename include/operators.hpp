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
  void setupOperators(int eType, int order, geo* inGeo, input* inParams);

  //! Setup operator for extrapolation from solution points to flux points
  void setupExtrapolatesptsfpts(vector<point> loc_spts, vector<point> loc_fpts, int eType, int order);

  //! Setup operator for calculation of gradient at the solution points
  void setupGradspts(vector<point> loc_spts, int eType, int order);

  //! Setup an interpolation operation between two sets of points using solution basis
  void setupInterpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp, int eType, int order);

  /*! Setup operator to calculate divergence of correction function at solution points
   *  based upon the normal flux correction at the flux points */
  void setup_correction(vector<point> loc_spts, vector<point> loc_fpts, int eType, int order);

  // Create a map<int,double*> (?) to get access to the correct operator
  // i.e. somthing like: div_flux_spts_tri = oper.get_oper_div[TRI]

  void apply_grad_spts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts);

  void apply_gradF_spts(vector<matrix<double>> &F_spts, vector<vector<matrix<double>>> &dF_spts);

  void apply_divF_spts(vector<matrix<double>> &F_spts, matrix<double> &divF_spts);

  void apply_spts_fpts(matrix<double> &U_spts, matrix<double> &U_fpts);

  const matrix<double>& get_oper_div_spts();
  const matrix<double>& get_oper_spts_fpts();

  //map<int,matrix<double>*> get_oper_div_spts;
  map<int,matrix<double>*> get_oper_grad_spts;
  //map<int,matrix<double>*> get_oper_spts_fpts;
  map<int,matrix<double>*> get_oper_correct;

private:
  geo *Geo;
  input *params;
  int nDims, eType, order;

  matrix<double> opp_spts_to_fpts;
  vector<matrix<double>> opp_grad_spts;
  matrix<double> opp_div_spts;
  matrix<double> opp_correction;
};
