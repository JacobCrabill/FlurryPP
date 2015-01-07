/*!
 * \file ele.cpp
 * \brief ele class definition
 *
 * Each elements stores its solution and basic properties, like
 * element type, vertices, and polynomial order
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill
 *
 */
#pragma once

#include "../include/ele.hpp"
#include "../include/polynomials.hpp"
#include "../include/flux.hpp"

using namespace std;

ele::ele(int in_eType, int in_order, int in_ID, vector<int> &in_nodes, mesh *in_Mesh)
{
  eType = in_eType;
  order = in_order;
  ID = in_ID;
  Mesh = in_Mesh;

  nodes.clear();
  for (auto &n: in_nodes)
    nodes.push_back(n);

}

ele::initialize()
{
  int fpt, face;

  switch(eType) {
  case(TRI):
    for (fpt=0; fpt<nFpts; fpt++) {
      face = fpt%(nFpts/3);
      switch(face) {
      case(0):
        tNorm_fpts[fpt][0] = 0;
        tNorm_fpts[fpt][1] = -1;
      case(1):
        tNorm_fpts[fpt][0] = 1.0/sqrt(2);
        tNorm_fpts[fpt][1] = 1.0/sqrt(2);
      case(2):
        tNorm_fpts[fpt][0] = -1;
        tNorm_fpts[fpt][1] = 0;
      }
    }
  case(QUAD):
    for (fpt=0; fpt<nFpts; fpt++) {
      face = fpt%(nFpts/3);
      switch(face) {
      case(0):
        tNorm_fpts[fpt][0] = 0;
        tNorm_fpts[fpt][1] = -1;
      case(1):
        tNorm_fpts[fpt][0] = 1;
        tNorm_fpts[fpt][1] = 0;
      case(2):
        tNorm_fpts[fpt][0] = 0;
        tNorm_fpts[fpt][1] = 1;
      case(2):
        tNorm_fpts[fpt][0] = -1;
        tNorm_fpts[fpt][1] = 0;
      }
    }
  }
}

void ele::calcTransforms(void)
{
  int spt, fpt, dim1, dim2;

  vector<vector<double>> dtemp(nDims,vector<double>(nDims));

  if (Jac_spts.size() != nSpts) Jac_spts.resize(nSpts);
  if (Jac_fpts.size() != nFpts) Jac_fpts.resize(nFpts);

  /* --- Calculate Transformation at Solution Points --- */
  for (spt=0; spt<nSpts; spt++) {
    // Calculate shape derivatives [in the future, should pre-calculate & store]
    switch(eType) {
    case TRI:
      dshape_tri(loc_spts[spt], dtemp);

    case QUAD:
      dshape_quad(loc_spts[spt], dtemp);

    default:
      FatalError("Element type not yet implemented.")
    }

    for (int i=0; i<nNodes; i++) {
      for (dim1=0; dim1<nDims; dim1++) {
        for (dim2=0; dim2<nDims; dim2++) {
          Jac_spts[spt][dim1][dim2] += dtemp[i][dim2]*nodes[i][dim1];
        }
      }
    }

    if (nDims==2) {
      detJac_spts[spt] = Jac_spts[spt][0][0]*Jac_spts[spt][1][1]-Jac_spts[spt][1][0]*Jac_spts[spt][0][1];
    }
    if (detJac_spts[spt]<0) FatalError("Negative Jacobian at solution points.");
  }

  /* --- Calculate Transformation at Flux Points --- */
  for (fpt=0; fpt<nFpts; fpt++) {
    // Calculate shape derivatives [in the future, should pre-calculate & store]
    switch(eType) {
    case TRI:
      dshape_tri(loc_fpts[fpt], dtemp);

    case QUAD:
      dshape_quad(loc_fpts[fpt], dtemp);

    default:
      FatalError("Element type not yet implemented.")
    }

    for (int i=0; i<nNodes; i++) {
      for (dim1=0; dim1<nDims; dim1++) {
        for (dim2=0; dim2<nDims; dim2++) {
          Jac_fpts[fpt][dim1][dim2] += dtemp[i][dim2]*nodes[i][dim1];
        }
      }
    }

    if (nDims==2) {
      detJac_fpts[fpt] = Jac_fpts[fpt][0][0]*Jac_fpts[fpt][1][1]-Jac_fpts[fpt][1][0]*Jac_fpts[fpt][0][1];
    }
    if (detJac_spts[spt]<0) FatalError("Negative Jacobian at solution points.");
  }
}

void ele::calcInviscidFlux_spts()
{
  for (int spt=0; spt<nSpts; spt++) {
    inviscidFlux(U_spts[spt], F_spts[spt], params);
  }
}

void ele::calcViscousFlux_spts()
{
  for (int spt=0; spt<nSpts; spt++) {
    viscousFlux(U_spts[spt], dU_spts[spt], F_spts[spt], params);
  }
}
