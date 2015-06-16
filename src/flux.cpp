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
}


void viscousFlux(double* U, matrix<double> &gradU, matrix<double> &Fvis, input *params)
{
  /* --- Calculate Primitives --- */
  double rho = U[0];
  double u   = U[1]/rho;
  double v   = U[2]/rho;
  double e   = U[3]/rho - 0.5*(u*u+v*v);

  /* --- Get Gradients --- */
  double dRho_dx	 = gradU(0,0);
  double dRhoU_dx = gradU(0,1);
  double dRhoV_dx = gradU(0,2);
  double dE_dx	   = gradU(0,3);

  double dRho_dy	 = gradU(1,0);
  double dRhoU_dy = gradU(1,1);
  double dRhoV_dy = gradU(1,2);
  double dE_dy	   = gradU(1,3);

  /* --- Calculate Viscosity --- */
  double rt_ratio = (params->gamma-1.0)*e/(params->rt_inf);
  double mu = (params->mu_inf)*pow(rt_ratio,1.5)*(1.+(params->c_sth))/(rt_ratio+(params->c_sth));
         mu+= params->fix_vis*(params->mu_inf - mu);

  double mu_t = 0;

  /* --- Calculate Gradients --- */
  double du_dx = (dRhoU_dx-dRho_dx*u)/rho;
  double du_dy = (dRhoU_dy-dRho_dy*u)/rho;

  double dv_dx = (dRhoV_dx-dRho_dx*v)/rho;
  double dv_dy = (dRhoV_dy-dRho_dy*v)/rho;

  double dK_dx = 0.5*(u*u+v*v)*dRho_dx+rho*(u*du_dx+v*dv_dx);
  double dK_dy = 0.5*(u*u+v*v)*dRho_dy+rho*(u*du_dy+v*dv_dy);

  double de_dx = (dE_dx-dK_dx-dRho_dx*e)/rho;
  double de_dy = (dE_dy-dK_dy-dRho_dy*e)/rho;

  double diag = (du_dx + dv_dy)/3.0;

  double tauxx = 2.0*(mu+mu_t)*(du_dx-diag);
  double tauxy = (mu+mu_t)*(du_dy + dv_dx);
  double tauyy = 2.0*(mu+mu_t)*(dv_dy-diag);

  /* --- Calculate Viscous Flux --- */
  Fvis(0,0) =  0.0;
  Fvis(0.1) = -tauxx;
  Fvis(0.2) = -tauxy;
  Fvis(0,3) = -(u*tauxx+v*tauxy+(mu/params->prandtl)*(params->gamma)*de_dx);

  Fvis(1,0) =  0.0;
  Fvis(1,1) = -tauxy;
  Fvis(1,2) = -tauyy;
  Fvis(1,3) = -(u*tauxy+v*tauyy+(mu/params->prandtl)*(params->gamma)*de_dy);
}

void viscousFluxAD(matrix<double> &gradU, matrix<double> &Fvis, input *params)
{
  Fvis(0,0) = params->diffD * gradU(0,0);
  Fvis(1,0) = params->diffD * gradU(1,0);
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
