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
    F[0][0] = params->advectVx*U[0];
    F[1][0] = params->advectVy*U[0];
  }
  else if (params->equation == NAVIER_STOKES) {
    double rho, u, v, p;
    rho = U[0];
    u = U[1]/rho;
    v = U[2]/rho;
    p = (params->gamma-1.0)*(U[3]-(0.5*rho*((u*u)+(v*v))));

    /* --- Assuming F has already been sized properly... --- */
    F[0][0] =  U[1];       F[1][0] =  U[2];
    F[0][1] =  U[1]*u+p;   F[1][1] =  U[1]*v;
    F[0][2] =  U[2]*u;     F[1][2] =  U[2]*v+p;
    F[0][3] = (U[3]+p)*u;  F[1][3] = (U[3]+p)*v;
  }
}


void viscousFlux(double* U, matrix<double> &gradU, matrix<double> &Fvis, input *params)
{
  double rho, u, v, E, e;

  double dRho_dx, dRhoU_dx, dRhoV_dx, dE_dx;
  double dRho_dy, dRhoU_dy, dRhoV_dy, dE_dy;

  double du_dx, du_dy, dv_dx, dv_dy, dK_dx, dK_dy, de_dx, de_dy;
  double diag, tauxx, tauxy, tauyy;
  double rt_ratio;

  double mu, mu_t, nu_tilde;
  double p,T,R;
  double inv_Re_c, Mach_c;
  double T_gas_non, S_gas_non;

  /* --- Calculate Primitives --- */
  rho = U[0];
  u   = U[1]/rho;
  v   = U[2]/rho;
  e   = U[3]/rho - 0.5*(u*u+v*v);

  /* --- Get Gradients --- */
  dRho_dx	 = gradU[0][0];
  dRhoU_dx = gradU[0][1];
  dRhoV_dx = gradU[0][2];
  dE_dx	   = gradU[0][3];

  dRho_dy	 = gradU[1][0];
  dRhoU_dy = gradU[1][1];
  dRhoV_dy = gradU[1][2];
  dE_dy	   = gradU[1][3];

  /* --- Calculate Viscosity --- */
  rt_ratio = (params->gamma-1.0)*e/(params->rt_inf);
  mu = (params->mu_inf)*pow(rt_ratio,1.5)*(1.+(params->c_sth))/(rt_ratio+(params->c_sth));
  mu = mu + params->fix_vis*(params->mu_inf - mu);

  /* --- Calculate Gradients --- */
  du_dx = (dRhoU_dx-dRho_dx*u)/rho;
  du_dy = (dRhoU_dy-dRho_dy*u)/rho;

  dv_dx = (dRhoV_dx-dRho_dx*v)/rho;
  dv_dy = (dRhoV_dy-dRho_dy*v)/rho;

  dK_dx = 0.5*(u*u+v*v)*dRho_dx+rho*(u*du_dx+v*dv_dx);
  dK_dy = 0.5*(u*u+v*v)*dRho_dy+rho*(u*du_dy+v*dv_dy);

  de_dx = (dE_dx-dK_dx-dRho_dx*e)/rho;
  de_dy = (dE_dy-dK_dy-dRho_dy*e)/rho;

  diag = (du_dx + dv_dy)/3.0;

  tauxx = 2.0*(mu+mu_t)*(du_dx-diag);
  tauxy = (mu+mu_t)*(du_dy + dv_dx);
  tauyy = 2.0*(mu+mu_t)*(dv_dy-diag);

  /* --- Calculate Viscous Flux --- */
  Fvis[0][0] =  0.0;
  Fvis[0][1] = -tauxx;
  Fvis[0][2] = -tauxy;
  Fvis[0][3] = -(u*tauxx+v*tauxy+(mu/params->prandtl)*(params->gamma)*de_dx);

  Fvis[1][0] =  0.0;
  Fvis[1][1] = -tauxy;
  Fvis[1][2] = -tauyy;
  Fvis[1][3] = -(u*tauxy+v*tauyy+(mu/params->prandtl)*(params->gamma)*de_dy);
}

//void rusanovFlux(vector<double> &UL, vector<double> &UR, vector<vector<double*>> &FL, vector<vector<double*>> &FR, vector<double> &norm, vector<double> &Fn, input *params)
void rusanovFlux(double* UL, double* UR, matrix<double> &FL, matrix<double> &FR, double* norm, double* Fn, input *params)
{
  int i, j;
  double rhoL, uL, vL, wL, pL, vnL=0.0;
  double rhoR, uR, vR, wR, pR, vnR=0.0;
  double csqL, csqR, eigL, eigR, eig;

  vector<double> FnL(params->nFields,0.0);
  vector<double> FnR(params->nFields,0.0);

  // Get primitive variables
  rhoL = UL[0];     rhoR = UR[0];
  uL = UL[1]/rhoL;  uR = UR[1]/rhoR;
  vL = UL[2]/rhoL;  vR = UR[2]/rhoR;

  // Calculate pressure
  if (params->nDims==2) {
    pL = (params->gamma-1.0)*(UL[3]-rhoL*(uL*uL+vL*vL));
    pR = (params->gamma-1.0)*(UR[3]-rhoR*(uR*uR+vR*vR));
  }else if (params->nDims==3) {
    wL = UL[3]/UL[0];   wR = UR[3]/UR[0];
    pL = (params->gamma-1.0)*(UL[3]-rhoL*(uL*uL+vL*vL+wL*wL));
    pR = (params->gamma-1.0)*(UR[3]-rhoR*(uR*uR+vR*vR+wR*wR));
  }

  // Get normal fluxes, normal velocities
  for (j=0; j<params->nDims; j++) {
    vnL += norm[j]*UL[j+1]/rhoL;
    vnR += norm[j]*UR[j+1]/rhoR;
    for (i=0; i<params->nFields; i++) {
      /*FnL[i] += norm[j]*(*FL[j][i]);
      FnR[i] += norm[j]*(*FR[j][i]);*/
      FnL[i] += norm[j]*FL[j][i];
      FnR[i] += norm[j]*FR[j][i];
    }
  }

  // Get maximum eigenvalue for diffusion coefficient
  csqL = max(params->gamma*pL/rhoL,0.0);
  csqR = max(params->gamma*pR/rhoR,0.0);
  eigL = fabs(vnL) + sqrt(csqL);
  eigR = fabs(vnR) + sqrt(csqR);
  eig = max(eigL,eigR);

  // Calculate Rusanov flux
  for (i=0; i<params->nFields; i++) {
    Fn[i] = 0.5*(FnL[i]+FnR[i] - eig*(UR[i]-UL[i]));
  }
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
      Fn[i] = 0.5*(FL[j][i]+FR[j][i])*norm[j];
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


void laxFriedrichsFlux(double* uL, double* uR, double* norm, double* Fn, input *params)
{
  if (params->equation == ADVECTION_DIFFUSION) {
    Fn[0] = params->advectVx*0.5*norm[0]*(uL[0]+uR[0])
          + params->advectVy*0.5*norm[1]*(uL[0]+uR[0]);

    double uAvg, uDiff, vNorm;

    uAvg = 0.5*(uL[0] + uR[0]);
    uDiff = uL[0] - uR[0];

    vNorm = params->advectVx*norm[0] + params->advectVy*norm[1];

    Fn[0] = vNorm*uAvg + 0.5*params->lambda*abs(vNorm)*uDiff;
  }
  else if (params->equation == NAVIER_STOKES) {
    FatalError("laxFlux not supported for Navier-Stokes simulations.");
  }
}

void ldgFlux(double* uL, double* uR, matrix<double> &gradU_L, matrix<double> &gradU_R, double* Fn, input *params)
{
  FatalError("LDG flux not implemented just yet.  Go to flux.cpp and do it!!");
}


void roeFlux(double *uL, double *uR, double* norm, double* Fn, input *params)
{
  double pL,pR;
  double hL, hR;
  double sq_rho,rrho,hm,usq,am,am_sq,unm,vgn;
  double lambda0,lambdaP,lambdaM;
  double rhoUnL, rhoUnR,eps;
  double a1,a2,a3,a4,a5,a6,aL1,bL1;

  double gamma = params->gamma;
  int nDims = params->nDims;

  array<double,3> vL, vR, um;
  array<double,5> du;

  // velocities
  for (int i=0;i<nDims;i++)  {
    vL[i] = uL[i+1]/uL[0];
    vR[i] = uR[i+1]/uR[0];
  }

  if (nDims==2) {
    pL = (gamma-1.0)*(uL[3] - (0.5*uL[0]*(vL[0]*vL[0]+vL[1]*vL[1])));
    pR = (gamma-1.0)*(uR[3] - (0.5*uR[0]*(vR[0]*vR[0]+vR[1]*vR[1])));
  }
  else
    FatalError("Roe not implemented in 3D");

  hL = (uL[nDims+1]+pL)/uL[0];
  hR = (uR[nDims+1]+pR)/uR[0];

  sq_rho = sqrt(uR[0]/uL[0]);
  rrho = 1./(sq_rho+1.);

  for (int i=0; i<nDims; i++)
    um[i] = rrho*(vL[i]+sq_rho*vR[i]);

  hm = rrho*(hL + sq_rho*hR);


  usq=0.;
  for (int i=0; i<nDims; i++)
    usq += 0.5*um[i]*um[i];

  am_sq = (gamma-1.)*(hm-usq);
  am = sqrt(am_sq);

  unm = 0.;
  vgn = 0.;
  for (int i=0;i<nDims;i++) {
    unm += um[i]*norm[i];
    //vgn += v_g(i)*norm[i];
  }

  // Compute Euler flux (first part)
  rhoUnL = 0.;
  rhoUnR = 0.;
  for (int i=0;i<nDims;i++)
  {
    rhoUnL += uL[i+1]*norm[i];
    rhoUnR += uR[i+1]*norm[i];
  }

  if (nDims==2)
  {
    Fn[0] = rhoUnL + rhoUnR;
    Fn[1] = rhoUnL*vL[0] + rhoUnR*vR[0] + (pL+pR)*norm[0];
    Fn[2] = rhoUnL*vL[1] + rhoUnR*vR[1] + (pL+pR)*norm[1];
    Fn[3] = rhoUnL*hL   +rhoUnR*hR;

  }
  else
    FatalError("Roe not implemented in 3D");

  for (int i=0;i<params->nFields;i++)
  {
    du[i] = uR[i]-uL[i];
  }

  lambda0 = abs(unm-vgn);
  lambdaP = abs(unm-vgn+am);
  lambdaM = abs(unm-vgn-am);

  // Entropy fix
  eps = 0.5*(abs(rhoUnL/uL[0]-rhoUnR/uR[0])+ abs(sqrt(gamma*pL/uL[0])-sqrt(gamma*pR/uR[0])));
  if(lambda0 < 2.*eps)
    lambda0 = 0.25*lambda0*lambda0/eps + eps;
  if(lambdaP < 2.*eps)
    lambdaP = 0.25*lambdaP*lambdaP/eps + eps;
  if(lambdaM < 2.*eps)
    lambdaM = 0.25*lambdaM*lambdaM/eps + eps;

  a2 = 0.5*(lambdaP+lambdaM)-lambda0;
  a3 = 0.5*(lambdaP-lambdaM)/am;
  a1 = a2*(gamma-1.)/am_sq;
  a4 = a3*(gamma-1.);

  if (nDims==2) {
    a5 = usq*du[0]-um[0]*du[1]-um[1]*du[2]+du[3];
    a6 = unm*du[0]-norm[0]*du[1]-norm[1]*du[2];
  }
  else if (nDims==3) {
    a5 = usq*du[0]-um[0]*du[1]-um[1]*du[2]-um[2]*du[3]+du[4];
    a6 = unm*du[0]-norm[0]*du[1]-norm[1]*du[2]-norm[2]*du[3];
  }

  aL1 = a1*a5 - a3*a6;
  bL1 = a4*a5 - a2*a6;

  // Compute Euler flux (second part)
  if (nDims==2) {
    Fn[0] -= lambda0*du[0]+aL1;
    Fn[1] -= lambda0*du[1]+aL1*um[0]+bL1*norm[0];
    Fn[2] -= lambda0*du[2]+aL1*um[1]+bL1*norm[1];
    Fn[3] -= lambda0*du[3]+aL1*hm   +bL1*unm;

  }
  else if (nDims==3) {
    Fn[0] -= lambda0*du[0]+aL1;
    Fn[1] -= lambda0*du[1]+aL1*um[0]+bL1*norm[0];
    Fn[2] -= lambda0*du[2]+aL1*um[1]+bL1*norm[1];
    Fn[3] -= lambda0*du[3]+aL1*um[2]+bL1*norm[2];
    Fn[4] -= lambda0*du[4]+aL1*hm   +bL1*unm;
  }

  for (int i=0;i<params->nFields;i++) {
    Fn[i] =  0.5*Fn[i];// - 0.5*vgn*(uR[i]+uL[i]);
  }

  //FatalError("Roe flux not implemented just yet.  Go to flux.cpp and do it!!");

}
