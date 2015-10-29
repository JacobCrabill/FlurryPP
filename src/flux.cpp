/*!
 * \file flux.cpp
 * \brief Functions for calculation of fluxes
 *
 * Includes functions for Euler and Navier-Stokes fluxes, as well as common flux
 * / Riemann solve routines
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "../include/flux.hpp"

#include <array>
#include <vector>

void inviscidFlux(double* U, matrix<double> &F, input *params)
{
  /* --- Note: Flux matrix expected to be <nDims x nFields> --- */
  if (params->equation == ADVECTION_DIFFUSION) {
    F(0,0) = params->advectVx*U[0];
    F(1,0) = params->advectVy*U[0];
    if (params->nDims == 3)
      F(2,0) = params->advectVz*U[0];
  }
  else if (params->equation == NAVIER_STOKES) {
    if (params->nDims == 2) {
      double rho = U[0];
      double u = U[1]/rho;
      double v = U[2]/rho;
      double p = (params->gamma-1.0)*(U[3]-(0.5*rho*((u*u)+(v*v))));

      /* --- Assuming F has already been sized properly... --- */
      F(0,0) =  U[1];       F(1,0) =  U[2];
      F(0,1) =  U[1]*u+p;   F(1,1) =  U[1]*v;
      F(0,2) =  U[2]*u;     F(1,2) =  U[2]*v+p;
      F(0,3) = (U[3]+p)*u;  F(1,3) = (U[3]+p)*v;
    }
    else if (params->nDims == 3) {
      double rho = U[0];
      double u = U[1]/rho;
      double v = U[2]/rho;
      double w = U[3]/rho;
      double p = (params->gamma-1.0)*(U[4]-(0.5*rho*((u*u)+(v*v)+(w*w))));

      /* --- Assuming F has already been sized properly... --- */
      F(0,0) =  U[1];       F(1,0) =  U[2];       F(2,0) =  U[3];
      F(0,1) =  U[1]*u+p;   F(1,1) =  U[1]*v;     F(2,1) =  U[1]*w;
      F(0,2) =  U[2]*u;     F(1,2) =  U[2]*v+p;   F(2,2) =  U[2]*w;
      F(0,3) =  U[3]*u;     F(1,3) =  U[3]*v;     F(2,3) =  U[3]*w+p;
      F(0,4) = (U[4]+p)*u;  F(1,4) = (U[4]+p)*v;  F(2,4) = (U[4]+p)*w;
    }
  }
}


void viscousFlux(double* U, matrix<double> &gradU, matrix<double> &Fvis, input *params)
{
  int nDims = params->nDims;

  /* --- Calculate Primitives --- */
  double rho = U[0];
  double u   = U[1]/rho;
  double v   = U[2]/rho;
  double e   = U[nDims+1]/rho - 0.5*(u*u+v*v);

  double w;
  if (nDims == 3) {
    w = U[3]/rho;
    e -= 0.5*(w*w);
  }
  else {
    w = 0.;
  }

  /* --- Get Gradients --- */
  double dRho_dx  = gradU(0,0);
  double dRhoU_dx = gradU(0,1);
  double dRhoV_dx = gradU(0,2);
  double dE_dx    = gradU(0,nDims+1);

  double dRho_dy  = gradU(1,0);
  double dRhoU_dy = gradU(1,1);
  double dRhoV_dy = gradU(1,2);
  double dE_dy	  = gradU(1,nDims+1);

  // 3D Derivatives
  double dRho_dz  = 0;
  double dRhoU_dz = 0;
  double dRhoV_dz = 0;
  double dRhoW_dx = 0;
  double dRhoW_dy = 0;
  double dRhoW_dz = 0;
  double dE_dz    = 0;
  if (nDims == 3) {
    dRho_dz	 = gradU(2,0);
    dRhoU_dz = gradU(2,1);
    dRhoV_dz = gradU(2,2);
    dRhoW_dx = gradU(0,3);
    dRhoW_dy = gradU(1,3);
    dRhoW_dz = gradU(2,3);
    dE_dz	   = gradU(2,4);
  }
//    /* --- Get Gradients --- */
//  double dRho_dx  = gradU[0];
//  double dRhoU_dx = gradU[1];
//  double dRhoV_dx = gradU[2];
//  double dE_dx    = gradU[nDims+1];
//
//  double dRho_dy  = gradU[nFields+0];
//  double dRhoU_dy = gradU[nFields+1];
//  double dRhoV_dy = gradU[nFields+2]
//  double dE_dy	  = gradU[nFields+nDims+1];
//
//  // 3D Derivatives
//  double dRho_dz  = 0;
//  double dRhoU_dz = 0;
//  double dRhoV_dz = 0;
//  double dRhoW_dx = 0;
//  double dRhoW_dy = 0;
//  double dRhoW_dz = 0;
//  double dE_dz    = 0;
//  if (nDims == 3) {
//    dRho_dz	 = gradU[2*nFields+0];
//    dRhoU_dz = gradU[2*nFields+1];
//    dRhoV_dz = gradU[2*nFields+2];
//    dRhoW_dx = gradU[3];
//    dRhoW_dy = gradU[nFields+3];
//    dRhoW_dz = gradU[2*nFields+3];
//    dE_dz	   = gradU[2*nFields+4];
//  }

  /* --- Calculate Viscosity --- */
  double mu = params->mu_inf;
  if (!params->fixVis) {
    // Use Sutherland's Law
    double rt_ratio = (params->gamma-1.0)*e/(params->rt_inf);
    mu *= pow(rt_ratio,1.5)*(1.+(params->c_sth))/(rt_ratio+(params->c_sth));
  }

  /* --- Calculate Gradients --- */
  double du_dx = (dRhoU_dx-dRho_dx*u)/rho;
  double du_dy = (dRhoU_dy-dRho_dy*u)/rho;

  double dv_dx = (dRhoV_dx-dRho_dx*v)/rho;
  double dv_dy = (dRhoV_dy-dRho_dy*v)/rho;

  // 3D Derivatives
  double du_dz=0, dv_dz=0;
  double dw_dx=0, dw_dy=0;
  double dw_dz = 0;
  if (nDims == 3) {
    du_dz = (dRhoU_dz-dRho_dz*u)/rho;
    dv_dz = (dRhoV_dz-dRho_dz*v)/rho;

    dw_dx = (dRhoW_dx-dRho_dx*w)/rho;
    dw_dy = (dRhoW_dy-dRho_dy*w)/rho;
    dw_dz = (dRhoW_dz-dRho_dz*w)/rho;
  }

  double dK_dx, dK_dy, dK_dz;
  if (nDims == 2) {
    dK_dx = 0.5*(u*u+v*v)*dRho_dx+rho*(u*du_dx+v*dv_dx);
    dK_dy = 0.5*(u*u+v*v)*dRho_dy+rho*(u*du_dy+v*dv_dy);
    dK_dz = 0;
  }
  else {
    dK_dx = 0.5*(u*u+v*v+w*w)*dRho_dx+rho*(u*du_dx+v*dv_dx+w*dw_dx);
    dK_dy = 0.5*(u*u+v*v+w*w)*dRho_dy+rho*(u*du_dy+v*dv_dy+w*dw_dy);
    dK_dz = 0.5*(u*u+v*v+w*w)*dRho_dz+rho*(u*du_dz+v*dv_dz+w*dw_dz);
  }

  double de_dx = (dE_dx-dK_dx-dRho_dx*e)/rho;
  double de_dy = (dE_dy-dK_dy-dRho_dy*e)/rho;
  double de_dz = 0;
  if (nDims == 3)
    de_dz = (dE_dz-dK_dz-dRho_dz*e)/rho;

  double diag = (du_dx + dv_dy + dw_dz)/3.0;

  double tauxx = 2.0*mu*(du_dx-diag);
  double tauyy = 2.0*mu*(dv_dy-diag);

  double tauxy = mu*(du_dy + dv_dx);

  double tauxz = 0;
  double tauyz = 0;
  double tauzz = 0;
  if (nDims == 3) {
    tauxz = mu*(du_dz + dv_dx);
    tauyz = mu*(du_dz + dv_dy);
    tauzz = 2.0*mu*(dw_dz-diag);
  }

  /* --- Calculate Viscous Flux --- */
  Fvis(0,0) =  0.0;
  Fvis(0,1) = -tauxx;
  Fvis(0,2) = -tauxy;
  Fvis(0,nDims+1) = -(u*tauxx+v*tauxy+w*tauxz+(mu/params->prandtl)*(params->gamma)*de_dx);

  Fvis(1,0) =  0.0;
  Fvis(1,1) = -tauxy;
  Fvis(1,2) = -tauyy;
  Fvis(1,nDims+1) = -(u*tauxy+v*tauyy+w*tauyz+(mu/params->prandtl)*(params->gamma)*de_dy);

  if (nDims == 3) {
    Fvis(0,3) = -tauzz;
    Fvis(1,3) = -tauyz;

    Fvis(2,0) =  0.0;
    Fvis(2,1) = -tauxz;
    Fvis(2,2) = -tauyz;
    Fvis(2,3) = -tauzz;
    Fvis(2,4) = -(u*tauxz+v*tauyz+w*tauzz+(mu/params->prandtl)*(params->gamma)*de_dz);
  }
}

matrix<double> viscousStressTensor(double* U, matrix<double> &gradU, input *params)
{
  int nDims = params->nDims;

  matrix<double> tau(nDims,nDims);

  /* --- Calculate Primitives --- */
  double rho = U[0];
  double u   = U[1]/rho;
  double v   = U[2]/rho;
  double e   = U[nDims+1]/rho - 0.5*(u*u+v*v);

  double w;
  if (nDims == 3) {
    w = U[3]/rho;
    e -= 0.5*(w*w);
  }
  else {
    w = 0.;
  }

  /* --- Get Gradients --- */
  double dRho_dx  = gradU(0,0);
  double dRhoU_dx = gradU(0,1);
  double dRhoV_dx = gradU(0,2);
  double dE_dx    = gradU(0,nDims+1);

  double dRho_dy  = gradU(1,0);
  double dRhoU_dy = gradU(1,1);
  double dRhoV_dy = gradU(1,2);
  double dE_dy	  = gradU(1,nDims+1);

  // 3D Derivatives
  double dRho_dz  = 0;
  double dRhoU_dz = 0;
  double dRhoV_dz = 0;
  double dRhoW_dx = 0;
  double dRhoW_dy = 0;
  double dRhoW_dz = 0;
  double dE_dz    = 0;
  if (nDims == 3) {
    dRho_dz	 = gradU(2,0);
    dRhoU_dz = gradU(2,1);
    dRhoV_dz = gradU(2,2);
    dRhoW_dx = gradU(0,3);
    dRhoW_dy = gradU(1,3);
    dRhoW_dz = gradU(2,3);
    dE_dz	   = gradU(2,4);
  }

  /* --- Calculate Viscosity --- */
  double mu = params->mu_inf;
  if (!params->fixVis) {
    // Use Sutherland's Law
    double rt_ratio = (params->gamma-1.0)*e/(params->rt_inf);
    mu *= pow(rt_ratio,1.5)*(1.+(params->c_sth))/(rt_ratio+(params->c_sth));
  }

  /* --- Calculate Gradients --- */
  double du_dx = (dRhoU_dx-dRho_dx*u)/rho;
  double du_dy = (dRhoU_dy-dRho_dy*u)/rho;

  double dv_dx = (dRhoV_dx-dRho_dx*v)/rho;
  double dv_dy = (dRhoV_dy-dRho_dy*v)/rho;

  // 3D Derivatives
  double du_dz=0, dv_dz=0;
  double dw_dx=0, dw_dy=0;
  double dw_dz = 0;
  if (nDims == 3) {
    du_dz = (dRhoU_dz-dRho_dz*u)/rho;
    dv_dz = (dRhoV_dz-dRho_dz*v)/rho;

    dw_dx = (dRhoW_dx-dRho_dx*w)/rho;
    dw_dy = (dRhoW_dy-dRho_dy*w)/rho;
    dw_dz = (dRhoW_dz-dRho_dz*w)/rho;
  }

  double dK_dx, dK_dy, dK_dz;
  if (nDims == 2) {
    dK_dx = 0.5*(u*u+v*v)*dRho_dx+rho*(u*du_dx+v*dv_dx);
    dK_dy = 0.5*(u*u+v*v)*dRho_dy+rho*(u*du_dy+v*dv_dy);
    dK_dz = 0;
  }
  else {
    dK_dx = 0.5*(u*u+v*v+w*w)*dRho_dx+rho*(u*du_dx+v*dv_dx+w*dw_dx);
    dK_dy = 0.5*(u*u+v*v+w*w)*dRho_dy+rho*(u*du_dy+v*dv_dy+w*dw_dy);
    dK_dz = 0.5*(u*u+v*v+w*w)*dRho_dz+rho*(u*du_dz+v*dv_dz+w*dw_dz);
  }

  double de_dx = (dE_dx-dK_dx-dRho_dx*e)/rho;
  double de_dy = (dE_dy-dK_dy-dRho_dy*e)/rho;
  double de_dz = 0;
  if (nDims == 3)
    de_dz = (dE_dz-dK_dz-dRho_dz*e)/rho;

  double diag = (du_dx + dv_dy + dw_dz)/3.0;

  for (int i=0; i<nDims; i++)
    tau(i,i) = diag;

  double tauxx = 2.0*mu*(du_dx-diag);
  double tauyy = 2.0*mu*(dv_dy-diag);

  tau(0,0) = tauxx;
  tau(1,1) = tauyy;

  double tauxy = mu*(du_dy + dv_dx);

  tau(0,1) = tauxy;
  tau(1,0) = tauxy;

  double tauxz = 0;
  double tauyz = 0;
  double tauzz = 0;
  if (nDims == 3) {
    tauxz = mu*(du_dz + dv_dx);
    tauyz = mu*(du_dz + dv_dy);
    tauzz = 2.0*mu*(dw_dz-diag);

    tau(0,2) = tauxz;
    tau(1,2) = tauyz;
    tau(2,2) = tauzz;

    tau(2,0) = tauxz;
    tau(2,1) = tauyz;
  }

  return tau;
}

void viscousFluxAD(matrix<double> &gradU, matrix<double> &Fvis, input *params)
{
  Fvis(0,0) = -params->diffD * gradU(0,0);
  Fvis(1,0) = -params->diffD * gradU(1,0);
  if (params->nDims == 3)
    Fvis(2,0) = -params->diffD * gradU(2,0);
}

void centralFlux(double* uL, double* uR, double* norm, double* Fn, input *params)
{
  if (params->equation == ADVECTION_DIFFUSION) {
    Fn[0] = params->advectVx*0.5*norm[0]*(uL[0]+uR[0])
          + params->advectVy*0.5*norm[1]*(uL[0]+uR[0]);
  }
  else if (params->equation == NAVIER_STOKES) {
    FatalError("centralFlux not supported for Navier-Stokes simulations.");
  }
}

void centralFlux(matrix<double> &FL, matrix<double> &FR, double* norm, double* Fn, input *params)
{
  for (int i=0; i<params->nFields; i++) {
    Fn[i] = 0;
    for (int j=0; j<params->nDims; j++) {
      Fn[i] += 0.5*(FL[j][i]+FR[j][i])*norm[j];
    }
  }
}

void upwindFlux(double* uL, double* uR, double* norm, double* Fn, input *params)
{
  if (params->equation == ADVECTION_DIFFUSION) {
    double FL, FR, alpha;
    alpha = params->advectVx*norm[0]+params->advectVy*norm[1];
    FL = alpha*uL[0];
    FR = alpha*uR[0];
    Fn[0] = 0.5*( FR+FL - alpha*(FR-FL) );
  }
}


void ldgFlux(double* , double* , matrix<double> &, matrix<double> &, double* , input *)
{
  FatalError("LDG flux not implemented just yet.  Go to flux.cpp and do it!!");
}
