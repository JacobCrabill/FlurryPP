/*!
 * \file operators.cpp
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Fux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include <cmath>
 
#include "../include/operators.hpp"
 
/*
 * example usage:
 * get_oper_grad[order] = &oper_grad;
 */

oper::setup_operators(int eType, int order, mesh *inMesh)
{
  // Get access to basic data
  Mesh = inMesh;

  nDims = Mesh->nDims;

  vector<point> loc_spts = Mesh->get_loc_spts(eType,order);
  vector<point> loc_fpts = Mesh->get_loc_fpts(eType,order);

  // Set up each operator
  setup_extrapolate_spts_fpts(loc_spts, loc_fpts, eType, order);

  setup_grad_spts(loc_spts, eType, order);
}

oper::setup_extrapolate_spts_fpts(vector<point> loc_spts, vector<point> loc_fpts, int eType, int order)
{
  int spt, fpt, nSpts, nFpts;
  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  opp_spts_to_fpts.setup(nSpts,nFpts);

  for (fpt=0; fpt<nFpts; fpt++) {
    for (spt=0; spt<nSpts; spt++) {
      switch(eType) {
      case(TRI):
        opp_spts_to_fpts.data[fpt][spt] = eval_dubiner_basis_2d(loc_fpts[fpt],spt,order);

      case(QUAD):
        // First, get the i an j ID of the spt
        ispt = mod(spt,nSpts/(order+1));
        jspt = floor(spt/(order+1));
        opp_spts_to_fpts.data[fpt][spt] = Lagrange(loc_spts,loc_fpts[fpt].x,ispt) * Lagrange(loc_spts,loc_fpts[fpt].y,jspt);

      default:
        FatalError("Element type not yet supported.");
      }
    }
  }
}

oper::setup_interpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp, int eType, int order)
{
  int ptA, ptB, nPtsFrom, nPtsTo;
  nPtsFrom = pts_from.size();
  nPtsTo = pts_to.size();

  opp_interp.setup(nSpts,nFpts);

  for (ptB=0; ptB<nPtsTo; ptB++) {
    for (ptA=0; ptA<nPtsFrom; ptA++) {
      switch(eType) {
      case(TRI):
        opp_interp[ptB][ptA] = eval_dubiner_basis_2d(pts_to[ptB],ptA,order);

      case(QUAD):
        // First, get the i and j ID of the pt [tensor-product element]
        iptA = mod(ptA,nPtsFrom/(order+1));
        jptA = floor(ptA/(order+1));
        opp_interp[ptB][ptA] = Lagrange(pts_from,pts_to[ptB].x,iptA) * Lagrange(pts_from,pts_to[ptB].y,jptA);

      default:
        FatalError("Element type not yet supported.");
      }
    }
  }
}


oper::setup_grad_spts(vector<point> loc_spts, int eType, int order)
{
  int spt, nSpts, ispt, jspt;
  nSpts = loc_spts.size();

  opp_grad_spts.setup(nSpts,nDims);

  switch(eType) {
  case(TRI):
    for (spt=0; spt<nSpts; spt++) {
      opp_grad_spts[0][spt] = eval_dr_dubiner_basis_2d(loc_spts[spt],spt,order);
      opp_grad_spts[1][spt] = eval_ds_dubiner_basis_2d(loc_spts[spt],spt,order);
    }

  case(QUAD):
    vector<double> loc_spts_1D = Mesh->get_loc_spts_1D(order);
    for (spt=0; spt<nSpts; spt++) {
      ispt = mod(ptA,nPtsFrom/(order+1));
      jspt = floor(ptA/(order+1));
      opp_grad_spts[0][spt] = dLagrange(loc_spts_1D,loc_spts_1D[spt].x,iptA) * Lagrange(loc_spts_1D,loc_spts_1D[spt].y,jptA);
      opp_grad_spts[1][spt] = Lagrange(loc_spts_1D,loc_spts_1D[spt].x,iptA) * dLagrange(loc_spts_1D,loc_spts_1D[spt].y,jptA);
    }

  default:
    FatalError("Element type not yet supported.");
  }
}


void oper::apply_grad_spts(matrix<double> &U_spts, vector<matrix<double>> &dU_spts)
{
  opp_grad_spts.timesMatrix(U_spts,dU_spts);
}

void oper::apply_spts_fpts(matrix<double> &U_spts, matrix<double> &U_fpts)
{
  opp_spts_to_fpts.timesMatrix(U_spts,U_fpts);
}


const matrix<double> &oper::get_oper_div_spts()
{
  return opp_div_spts;
}
