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

#include "../include/ele.hpp"
#include "../include/polynomials.hpp"
#include "../include/flux.hpp"

using namespace std;

ele::ele()
{

}

ele::ele(int in_eType, int in_order, int in_ID, vector<point> &in_nodes, geo *in_Geo)
{
  eType = in_eType;
  order = in_order;
  ID = in_ID;
  Geo = in_Geo;

  nodes.clear();
  for (auto &n: in_nodes)
    nodes.push_back(n);

}

void ele::initialize(void)
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
        break;
      case(1):
        tNorm_fpts[fpt][0] = 1.0/sqrt(2);
        tNorm_fpts[fpt][1] = 1.0/sqrt(2);
        break;
      case(2):
        tNorm_fpts[fpt][0] = -1;
        tNorm_fpts[fpt][1] = 0;
        break;
      }
    }
    break;
  case(QUAD):
    for (fpt=0; fpt<nFpts; fpt++) {
      face = fpt%(nFpts/3);
      switch(face) {
      case(0):
        tNorm_fpts[fpt][0] = 0;
        tNorm_fpts[fpt][1] = -1;
        break;
      case(1):
        tNorm_fpts[fpt][0] = 1;
        tNorm_fpts[fpt][1] = 0;
        break;
      case(2):
        tNorm_fpts[fpt][0] = 0;
        tNorm_fpts[fpt][1] = 1;
        break;
      case(3):
        tNorm_fpts[fpt][0] = -1;
        tNorm_fpts[fpt][1] = 0;
        break;
      }
    }
    break;
  }
}

void ele::setup(input *inParams, geo *inGeo)
{
  /* --- Basic Stuff --- */
  params = inParams;
  Geo = inGeo;

  order = params->order;
  nDims = params->nDims;

  loc_spts = Geo->getLocSpts(eType,order);
  loc_fpts = Geo->getLocFpts(eType,order);

  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  pos_spts.resize(nSpts);

  if (params->equation == ADVECTION_DIFFUSION) {
    nFields = 1;
  }else if (params->equation == NAVIER_STOKES) {
    nFields = nDims + 2;
  }

  /* --- Setup all data arrays --- */
  U_spts.setup(nSpts,nFields);
  U_fpts.setup(nFpts,nFields);
  Fn_fpts.setup(nFpts,nFields);
  divF_spts.setup(nSpts,nFields);

  F_spts.resize(nSpts);
  F_fpts.resize(nFpts);
  dU_spts.resize(nSpts);
  dU_fpts.resize(nFpts);
  for (int spt=0; spt<nSpts; spt++) { // if I fix nDims (change to nSpts) code fails at line 22 of flurry.cpp!
    F_spts[spt].setup(nDims,nFields);
    dU_spts[spt].setup(nDims,nFields);
  }
  for (int fpt=0; fpt<nFpts; fpt++) {
    F_fpts[fpt].setup(nDims,nFields);
    dU_fpts[fpt].setup(nDims,nFields);
  }

  //dF_spts.setup(nDims,nDims);
  dF_spts.resize(nDims);
  for (int i=0; i<nDims; i++) {
    dF_spts[i].resize(nDims);
    for (int j=0; j<nDims; j++) {
      dF_spts[i][j].setup(nSpts,nFields);
    }
  }

  detJac_spts.resize(nSpts);
  detJac_fpts.resize(nFpts);
  Jac_spts.resize(nSpts);
  Jac_fpts.resize(nFpts);
  for (auto& spt:Jac_spts) spt.setup(nDims,nDims);
  for (auto& fpt:Jac_fpts) fpt.setup(nDims,nDims);

  norm_fpts.setup(nFpts,nDims);
  tNorm_fpts.setup(nFpts,nDims);

  /* --- Final Step: calculate physical->reference transforms --- */
  calcTransforms();

  calcPosSpts();
}

void ele::calcTransforms(void)
{
  matrix<double> dtemp(nDims,nDims);

  if (Jac_spts.size() != (uint)nSpts) Jac_spts.resize(nSpts);
  if (Jac_fpts.size() != (uint)nFpts) Jac_fpts.resize(nFpts);

  /* --- Calculate Transformation at Solution Points --- */
  for (int spt=0; spt<nSpts; spt++) {
    // Calculate shape derivatives [in the future, should pre-calculate & store]
    switch(eType) {
      case TRI:
        dshape_tri(loc_spts[spt], dtemp);
        break;
      case QUAD:
        dshape_quad(loc_spts[spt], dtemp);
        break;
      default:
        FatalError("Element type not yet implemented.")
    }

    for (int i=0; i<nNodes; i++) {
      for (int dim1=0; dim1<nDims; dim1++) {
        for (int dim2=0; dim2<nDims; dim2++) {
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
  for (int fpt=0; fpt<nFpts; fpt++) {
    // Calculate shape derivatives [in the future, should pre-calculate & store]
    switch(eType) {
      case TRI:
        dshape_tri(loc_fpts[fpt], dtemp);
        break;
      case QUAD:
        dshape_quad(loc_fpts[fpt], dtemp);
        break;
      default:
        FatalError("Element type not yet implemented.")
    }

    for (int i=0; i<nNodes; i++) {
      for (int dim1=0; dim1<nDims; dim1++) {
        for (int dim2=0; dim2<nDims; dim2++) {
          Jac_fpts[fpt][dim1][dim2] += dtemp[i][dim2]*nodes[i][dim1];
        }
      }
    }

    if (nDims==2) {
      detJac_fpts[fpt] = Jac_fpts[fpt][0][0]*Jac_fpts[fpt][1][1]-Jac_fpts[fpt][1][0]*Jac_fpts[fpt][0][1];
    }
    if (detJac_fpts[fpt]<0) FatalError("Negative Jacobian at solution points.");
  }
}

void ele::calcPosSpts(void)
{
  vector<double> shape;

  for (int spt=0; spt<nSpts; spt++) {
    getShape(spt, shape);
    pos_spts[spt].zero();
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        pos_spts[spt][dim] += shape[iv]*nodes[iv][dim];
      }
    }
  }
}

void ele::setInitialCondition()
{
  if (params->equation == NAVIER_STOKES) {
    double rho, vx, vy, p;
    double gamma = params->gamma;

    if (params->ic_type == 0) {
      /* --- Uniform "Freestream" solution --- */
      rho = params->rhoIC;
      vx = params->vxIC;
      vy = params->vyIC;
      p = params->pIC;
      for (int spt=0; spt<nSpts; spt++) {
        U_spts[spt][0] = rho;
        U_spts[spt][1] = rho * vx;
        U_spts[spt][2] = rho * vy;
        U_spts[spt][3] = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
    else if (params->ic_type == 1) {
      /* --- Isentropic Vortex of strength eps centered at (0,0) --- */
      double eps, f, x, y;
      for (int spt=0; spt<nSpts; spt++) {
        eps = 5.0;
        x = pos_spts[spt][0];
        y = pos_spts[spt][1];

        f = 1.0 - (x*x + y*y);

        rho = pow(1. - eps*eps*(gamma-1.)/(8.*gamma*pi*pi)*exp(f), 1.0/(gamma-1.0));
        vx = 1. - eps*y / (2.*pi) * exp(f/2.);
        vy = 1. + eps*x / (2.*pi) * exp(f/2.);
        p = pow(rho,gamma);

        U_spts[spt][0] = rho;
        U_spts[spt][1] = rho * vx;
        U_spts[spt][2] = rho * vy;
        U_spts[spt][3] = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    /* --- Simple Gaussian bump centered at (0,0) --- */
    double r2;
    for (int spt=0; spt<nSpts; spt++) {
      r2 = pos_spts[spt][0]*pos_spts[spt][0] + pos_spts[spt][1]*pos_spts[spt][1];
      U_spts[spt][0] = exp(-r2);
    }
  }
}

void ele::getShape(int spt, vector<double> &shape)
{
  if (eType == TRI) {
    shape_tri(loc_spts[spt], shape);
  }
  else if (eType == QUAD) {
    shape_quad(loc_spts[spt], shape);
  }
  else {
    FatalError("Element Type Not Supported.");
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

vector<double> ele::getPrimitives(uint spt)
{
  vector<double> V(nFields);

  if (params->equation == ADVECTION_DIFFUSION) {
    V[0] = U_spts[spt][0];
  }
  else if (params->equation == NAVIER_STOKES) {
    V[0] = U_spts[spt][0];
    V[1] = U_spts[spt][1]/V[0];
    V[2] = U_spts[spt][2]/V[0];
    // p = (gamma-1)*(rho*e - 1/2*rho*(magV^2))
    V[3] = (params->gamma-1)*(U_spts[spt][3] - (0.5*(U_spts[spt][1]*U_spts[spt][1] + U_spts[spt][2]*U_spts[spt][2])/U_spts[spt][0]));
  }

  return V;
}

vector<double> ele::getResidual(int normType)
{
  vector<double> res(nFields,0);

  for (int spt=0; spt<nSpts; spt++) {
    for (int i=0; i<nFields; i++) {
      if (normType == 1) {
        res[i] += abs(U_spts[spt][i]);
      }
      else if (normType == 2) {
        res[i] += U_spts[spt][i]*U_spts[spt][i];
      }
      else if (normType == 3) {
        // Infinity norm
        res[i] = max(abs(U_spts[spt][i]),res[i]);
      }
    }
  }

  return res;
}

point ele::getPosSpt(uint spt)
{
  return pos_spts[spt];
}
int ele::getNDims() const
{
  return nDims;
}

void ele::setNDims(int value)
{
  nDims = value;
}
int ele::getNFields() const
{
  return nFields;
}

void ele::setNFields(int value)
{
  nFields = value;
}
int ele::getNSpts() const
{
  return nSpts;
}

void ele::setNSpts(int value)
{
  nSpts = value;
}
int ele::getNFpts() const
{
  return nFpts;
}

void ele::setNFpts(int value)
{
  nFpts = value;
}




