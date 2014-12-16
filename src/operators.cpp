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

oper::setup_operators(int eleType, int order, mesh &_Mesh)
{
  vector<double> loc_spts = Mesh->get_loc_spts(eleType,order);
  vector<double> loc_fpts = Mesh->get_loc_fpts(eleType,order);
}

oper::setup_extrapolate_spts_fpts(vector<point> loc_spts, vector<point> loc_fpts, int order, int eType)
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

oper::setup_grad_spts(vector<point> loc_spts, int order, int eType)
{

}
