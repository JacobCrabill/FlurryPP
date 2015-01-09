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

void oper::setup_operators(int eType, int order, geo *inGeo)
{
  // Get access to basic data
  Geo = inGeo;

  nDims = Geo->nDims;

  vector<point> loc_spts = Geo->getLocSpts(eType,order);
  vector<point> loc_fpts = Geo->getLocFpts(eType,order);

  // Set up each operator
  setup_extrapolate_spts_fpts(loc_spts, loc_fpts, eType, order);

  setup_grad_spts(loc_spts, eType, order);
}

void oper::setup_extrapolate_spts_fpts(vector<point> loc_spts, vector<point> loc_fpts, int eType, int order)
{
  unsigned int spt, fpt, nSpts, nFpts, ispt, jspt;
  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  opp_spts_to_fpts.setup(nSpts,nFpts);

  for (fpt=0; fpt<nFpts; fpt++) {
    for (spt=0; spt<nSpts; spt++) {
      switch(eType) {
      case(TRI):
        opp_spts_to_fpts[fpt][spt] = eval_dubiner_basis_2d(loc_fpts[fpt],spt,order);
        break;
      case(QUAD): {
        vector<double> locSpts1D = Geo->getLocSpts1D(eType,order);
        // First, get the i an j ID of the spt
        ispt = spt%(nSpts/(order+1));
        jspt = floor(spt/(order+1));
        opp_spts_to_fpts[fpt][spt] = Lagrange(locSpts1D,loc_fpts[fpt].x,ispt) * Lagrange(locSpts1D,loc_fpts[fpt].y,jspt);
        break;
      }
      default:
        FatalError("Element type not yet supported.");
      }
    }
  }
}

void oper::setup_interpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp, int eType, int order)
{
  unsigned int ptA, ptB, nPtsFrom, nPtsTo, pt, iptA, jptA;
  nPtsFrom = pts_from.size();
  nPtsTo = pts_to.size();

  opp_interp.setup(nPtsFrom,nPtsTo);

  // Get 1D locations of points for arbitrary, potentially anisotropic tensor-product elements
  vector<double> locPts1Dx, locPts1Dy;
  if (eType==QUAD) {
    locPts1Dx.resize(order+1);
    locPts1Dy.resize(order+1);
    for (pt=0; pt<(unsigned int)order+1; pt++) {
      locPts1Dx[pt] = pts_from[pt].x;
      locPts1Dy[pt] = pts_from[pt*(order+1)].y;
    }
  }

  for (ptB=0; ptB<nPtsTo; ptB++) {
    for (ptA=0; ptA<nPtsFrom; ptA++) {
      switch(eType) {
      case(TRI):
        opp_interp[ptB][ptA] = eval_dubiner_basis_2d(pts_to[ptB],ptA,order);
        break;
      case(QUAD):
        // First, get the i and j ID of the pt [tensor-product element]
        iptA = ptA%(nPtsFrom/(order+1));
        jptA = floor(ptA/(order+1));
        opp_interp[ptB][ptA] = Lagrange(locPts1Dx,pts_to[ptB].x,iptA) * Lagrange(locPts1Dy,pts_to[ptB].y,jptA);
        break;
      default:
        FatalError("Element type not yet supported.");
      }
    }
  }
}


void oper::setup_grad_spts(vector<point> loc_spts, int eType, int order)
{
  unsigned int spt, nSpts, ispt, jspt;
  nSpts = loc_spts.size();

  opp_grad_spts.setup(nSpts,nDims);

  switch(eType) {
  case(TRI):
    for (spt=0; spt<nSpts; spt++) {
      opp_grad_spts[0][spt] = eval_dr_dubiner_basis_2d(loc_spts[spt],spt,order);
      opp_grad_spts[1][spt] = eval_ds_dubiner_basis_2d(loc_spts[spt],spt,order);
    }
    break;
  case(QUAD): {
    vector<double> loc_spts_1D = Geo->getLocSpts1D(eType,order);
    for (spt=0; spt<nSpts; spt++) {
      ispt = spt%(nSpts/(order+1));
      jspt = floor(spt/(order+1));
      opp_grad_spts[0][spt] = dLagrange(loc_spts_1D,loc_spts_1D[spt],ispt) * Lagrange(loc_spts_1D,loc_spts_1D[spt],jspt);
      opp_grad_spts[1][spt] = Lagrange(loc_spts_1D,loc_spts_1D[spt],ispt) * dLagrange(loc_spts_1D,loc_spts_1D[spt],jspt);
    }
    break;
  }
  default:
    FatalError("Element type not yet supported.");
  }
}


void oper::apply_grad_spts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts)
{
  for (unsigned int dim=0; dim<dU_spts.size(); dim++)
    opp_grad_spts.timesMatrix(U_spts,dU_spts[dim]);
}

void oper::apply_spts_fpts(matrix<double> &U_spts, matrix<double> &U_fpts)
{
  opp_spts_to_fpts.timesMatrix(U_spts,U_fpts);
}


const matrix<double> &oper::get_oper_div_spts()
{
  return opp_div_spts;
}
