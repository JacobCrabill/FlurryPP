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
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
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

  /* --- Calculate Viscosity --- */
  double rt_ratio = (params->gamma-1.0)*e/(params->rt_inf);
  double mu = (params->mu_inf)*pow(rt_ratio,1.5)*(1.+(params->c_sth))/(rt_ratio+(params->c_sth));
         mu+= params->fix_vis*(params->mu_inf - mu);

  double mu_t = 0.;

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

  double tauxx = 2.0*(mu+mu_t)*(du_dx-diag);
  double tauyy = 2.0*(mu+mu_t)*(dv_dy-diag);

  double tauxy = (mu+mu_t)*(du_dy + dv_dx);

  double tauxz = 0;
  double tauyz = 0;
  double tauzz = 0;
  if (nDims == 3) {
    tauxz = (mu+mu_t)*(du_dz + dv_dx);
    tauyz = (mu+mu_t)*(du_dz + dv_dy);
    tauzz = 2.0*(mu+mu_t)*(dw_dz-diag);
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
