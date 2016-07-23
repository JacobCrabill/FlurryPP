/*!
 * \file face.cpp
 * \brief Abstract parent class to handle interface flux calculations at all faces
 *
 * All type-specific face classes will be derived from this class
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

#include "face.hpp"

#include "flux.hpp"
#include "ele.hpp"

void face::initialize(shared_ptr<ele> &eL, shared_ptr<ele> &eR, int gID, int locF_L, faceInfo myInfo, input *params)
{
  ID = gID;

  this->eL = eL;
  this->eR = eR;
  this->locF_L = locF_L;
  this->myInfo = myInfo;
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;

  // Setup temporary vectors for later use
  //tempFL.setup(nDims,nFields);
  //tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);
}

void face::setupFace(void)
{
  if (nDims == 2)
    nFptsL = eL->order+1;
  else
    nFptsL = (eL->order+1)*(eL->order+1);

  // Find first/last ID of fpts on given face
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.setup(nFptsL,nFields);
  UR.setup(nFptsL,nFields);
  FnL.resize(nFptsL);
  Fn.setup(nFptsL,nFields);
  normL.setup(nFptsL,nDims);
  dAL.resize(nFptsL);
  //detJacL.resize(nFptsL);
  waveSp.resize(nFptsL);

  Fn.initializeToZero();

  if (params->viscous) {
    UC.setup(nFptsL,nFields);
    dUcL.setup(nFptsL,nFields);
    // just a placeholder.  Need to properly size/reorder dimensions later.
    gradUL.resize(nFptsL);
    gradUR.resize(nFptsL);
    for (auto &dU:gradUL) dU.setup(nDims,nFields);
    for (auto &dU:gradUR) dU.setup(nDims,nFields);
  }

  if (params->motion) {
    Vg.setup(nFptsL,nDims);
  }

  getPointers();

  this->setupRightState();

  this->getPointersRight();
}

void face::getPointers(void)
{
  // Get access to data at left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    FnL[fpt] = &(eL->Fn_fpts(i,0));
    waveSp[fpt] = &(eL->waveSp_fpts[i]);

    if (params->viscous) {
      for  (int k = 0; k < nFields; k++)
        dUcL(fpt,k) = &(eL->dUc_fpts(i,k));
    }

    fpt++;
  }
}

void face::getLeftState()
{
  // Get data from left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    for (int j=0; j<nFields; j++) {
      UL(fpt,j) = (eL->U_fpts(i,j));
    }

    /* For dynamic grids (besides rigid translation), need to update
     * geometry-related data on every iteration, not just during setup */
    if (isNew || (params->motion != 0 && params->motion != 4)) {
      for (int dim=0; dim<nDims; dim++) {
        normL(fpt,dim) = (eL->norm_fpts(i,dim));
      }
      dAL[fpt] = (eL->dA_fpts(i));
    }

    if (params->motion) {
      for (int dim=0; dim<nDims; dim++)
        Vg(fpt,dim) = eL->gridVel_fpts(i,dim);
    }

    fpt++;
  }

  isNew = false;
}

void face::getLeftGradient()
{
  // Get data from left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    if (params->viscous) {
      for (int dim=0; dim<nDims; dim++)
        for (int j=0; j<nFields; j++)
          gradUL[fpt](dim,j) = eL->dU_fpts(dim,i,j);
    }

    fpt++;
  }
}

void face::calcInviscidFlux(void)
{
  if (!isMPI)
    getLeftState();
  this->getRightState(); // <-- makes this more general for all face types, and allows face memory to be contiguous

  // Calculate common inviscid flux at flux points
  if (params->equation == ADVECTION_DIFFUSION) {
    laxFriedrichsFlux();
  }
  else if (params->equation == NAVIER_STOKES) {
    if (isBnd) {
      //centralFluxBound();
      rusanovFlux();
    } else {
      if (params->riemannType==0) {
        rusanovFlux();
      }
      else if (params->riemannType==1) {
        roeFlux();
      }
    }
  }

  // Transform normal flux using edge Jacobian and put into ele's memory
  for (int i=0; i<nFptsL; i++)
    for (int j=0; j<nFields; j++)
      FnL[i][j] =  Fn(i,j)*dAL[i];

  if (params->viscous) {

    ldgSolution();

    // Still assuming nFptsL = nFptsR
    for (int i=0; i<nFptsL; i++) {
      for (int j=0; j<nFields; j++) {
        *dUcL(i,j) = UC(i,j) - UL(i,j);
      }
    }
  }

  this->setRightStateFlux();

  if (params->viscous)
    this->setRightStateSolution();
}

void face::calcViscousFlux(void)
{
  if (!isMPI)
    getLeftGradient();
  this->getRightGradient();

  if (params->equation == NAVIER_STOKES) {
    for (int fpt=0; fpt<nFptsL; fpt++) {
      // Calculte common viscous flux at flux points [LDG numerical flux]
      double Fc[nDims][nFields];

      if (isBnd) {
        if (isBnd > 1) {
          // Adiabatic wall boundary condition (Neumann-type BC)
          // Calculate discontinuous viscous flux at right flux points
          viscousFlux(UR[fpt], gradUR[fpt], tempFR, params);
          for (int dim=0; dim<nDims; dim++) {
            for (int k=0; k<nFields; k++) {
              Fc[dim][k] = tempFR[dim][k] + params->tau*normL(fpt,dim)*(UL(fpt,k) - UR(fpt,k));
            }
          }
        }
        else {
          // All other boundary conditions (Dirichlet-type BC)
          // Calculate discontinuous viscous flux at left flux points
          viscousFlux(UL[fpt], gradUL[fpt], tempFL, params);
          for (int dim=0; dim<nDims; dim++) {
            for (int k=0; k<nFields; k++) {
              Fc[dim][k] = tempFL[dim][k] + params->tau*normL(fpt,dim)*(UL(fpt,k) - UR(fpt,k));
            }
          }
        }
      }
      else {
        // All general interior-type faces (interior, MPI, overset)
        viscousFlux(UL[fpt], gradUL[fpt], tempFL, params);
        viscousFlux(UR[fpt], gradUR[fpt], tempFR, params);
        double penFact = params->penFact;
        if (nDims == 2) {
          if ( normL(fpt,0)+normL(fpt,1) < 0 )
            penFact = -params->penFact;
        }
        else if (nDims == 3) {
          if (normL(fpt,0)+normL(fpt,1)+sqrt(2.)*normL(fpt,2) < 0) {
            penFact = -params->penFact;
          }
        }

        double normX = normL(fpt,0);
        double normY = normL(fpt,1);
        double normZ = 0;
        if (nDims == 3)
          normZ = normL(fpt,2);

        if (nDims == 2) {
          for(int k=0; k<nFields; k++) {
            double dF0 = tempFL[0][k] - tempFR[0][k];
            double dF1 = tempFL[1][k] - tempFR[1][k];
            double dU = UL(fpt,k) - UR(fpt,k);
            Fc[0][k] = 0.5*(tempFL[0][k] + tempFR[0][k]) + penFact*normX*( normX*(dF0) + normY*(dF1) ) + params->tau*normX*(dU);
            Fc[1][k] = 0.5*(tempFL[1][k] + tempFR[1][k]) + penFact*normY*( normX*(dF0) + normY*(dF1) ) + params->tau*normY*(dU);
          }
        }
        else if (nDims == 3) {
          for(int k=0; k<nFields; k++) {
            double dF0 = tempFL[0][k] - tempFR[0][k];
            double dF1 = tempFL[1][k] - tempFR[1][k];
            double dF2 = tempFL[1][k] - tempFR[2][k];
            double dU = UL(fpt,k) - UR(fpt,k);
            Fc[0][k] = 0.5*(tempFL[0][k] + tempFR[0][k]) + penFact*normX*( normX*(dF0) + normY*(dF1) + normZ*(dF2) ) + params->tau*normX*(dU);
            Fc[1][k] = 0.5*(tempFL[1][k] + tempFR[1][k]) + penFact*normY*( normX*(dF0) + normY*(dF1) + normZ*(dF2) ) + params->tau*normY*(dU);
            Fc[2][k] = 0.5*(tempFL[2][k] + tempFR[2][k]) + penFact*normZ*( normX*(dF0) + normY*(dF1) + normZ*(dF2) ) + params->tau*normZ*(dU);
          }
        }
      }

      // calculate normal flux from discontinuous solution at flux points
      for (int dim=0; dim<nDims; dim++)
        for(int k=0; k<nFields; k++)
          Fn(fpt,k) += Fc[dim][k]*normL(fpt,dim);
    }
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    for (int fpt=0; fpt<nFptsL; fpt++) {
      viscousFluxAD(gradUL[fpt], tempFL, params);
      viscousFluxAD(gradUR[fpt], tempFR, params);

      double tempFn[1];
      centralFlux(tempFL, tempFR, normL[fpt], tempFn, params);

      Fn(fpt,0) += tempFn[0];
    }
  }

  // Transform normal flux using edge Jacobian and put into ele's memory
  for (int i=0; i<nFptsL; i++) {
    for (int j=0; j<nFields; j++)
      FnL[i][j] =  Fn(i,j)*dAL[i];
  }

  this->setRightStateFlux();
}


void face::rusanovFlux(void)
{
  for (int fpt=0; fpt<nFptsL; fpt++) {
    inviscidFlux(UL[fpt],tempFL,params);
    inviscidFlux(UR[fpt],tempFR,params);

    double tempFnL[5] = {0,0,0,0,0};
    double tempFnR[5] = {0,0,0,0,0};

    // Get primitive variables
    double rhoL = UL(fpt,0);     double rhoR = UR(fpt,0);
    double uL = UL(fpt,1)/rhoL;  double uR = UR(fpt,1)/rhoR;
    double vL = UL(fpt,2)/rhoL;  double vR = UR(fpt,2)/rhoR;

    double wL, pL, vnL=0.;
    double wR, pR, vnR=0.;
    double vgn=0.;

    // Calculate pressure
    if (params->nDims==2) {
      pL = (params->gamma-1.0)*(UL(fpt,3)-0.5*rhoL*(uL*uL+vL*vL));
      pR = (params->gamma-1.0)*(UR(fpt,3)-0.5*rhoR*(uR*uR+vR*vR));
    }
    else {
      wL = UL(fpt,3)/rhoL;   wR = UR(fpt,3)/rhoR;
      pL = (params->gamma-1.0)*(UL(fpt,4)-0.5*rhoL*(uL*uL+vL*vL+wL*wL));
      pR = (params->gamma-1.0)*(UR(fpt,4)-0.5*rhoR*(uR*uR+vR*vR+wR*wR));
    }

    // Get normal fluxes, normal velocities
    for (int dim=0; dim<params->nDims; dim++) {
      vnL += normL(fpt,dim)*UL(fpt,dim+1)/rhoL;
      vnR += normL(fpt,dim)*UR(fpt,dim+1)/rhoR;
      if (params->motion)
        vgn += normL(fpt,dim)*Vg(fpt,dim);
      for (int i=0; i<params->nFields; i++) {
        tempFnL[i] += normL(fpt,dim)*tempFL[dim][i];
        tempFnR[i] += normL(fpt,dim)*tempFR[dim][i];
      }
    }

    // Get maximum eigenvalue for diffusion coefficient
    double csqL = max(params->gamma*pL/rhoL,0.0);
    double csqR = max(params->gamma*pR/rhoR,0.0);
    double eigL = std::fabs(vnL) + sqrt(csqL);
    double eigR = std::fabs(vnR) + sqrt(csqR);
    double eig  = max(eigL,eigR);

    // Calculate Rusanov flux
    for (int i=0; i<params->nFields; i++) {
      Fn(fpt,i) = 0.5*(tempFnL[i]+tempFnR[i] - eig*(UR(fpt,i)-UL(fpt,i)));
    }

    // Store wave speed for calculation of allowable dt
    if (params->motion) {
      eigL = std::fabs(vnL-vgn) + sqrt(csqL);
      eigR = std::fabs(vnR-vgn) + sqrt(csqR);
    }
    *waveSp[fpt] = max(eigL,eigR);
  }
}

void face::roeFlux(void)
{
  double gamma = params->gamma;
  int nDims = params->nDims;

  if (nDims == 3) FatalError("Roe not implemented in 3D");

  for (int fpt=0; fpt<nFptsL; fpt++) {
    double vL[3], vR[3], um[3];
    double du[5];

    // velocities
    for (int i=0;i<nDims;i++)  {
      vL[i] = UL(fpt,i+1)/UL(fpt,0);
      vR[i] = UR(fpt,i+1)/UR(fpt,0);
    }

    double pL, pR;
    if (nDims == 2) {
      pL = (gamma-1.0)*(UL(fpt,3) - (0.5*UL(fpt,0)*(vL[0]*vL[0]+vL[1]*vL[1])));
      pR = (gamma-1.0)*(UR(fpt,3) - (0.5*UR(fpt,0)*(vR[0]*vR[0]+vR[1]*vR[1])));
    }
    else {
      pL = (gamma-1.0)*(UL(fpt,4) - (0.5*UL(fpt,0)*(vL[0]*vL[0]+vL[1]*vL[1]+vL[2]*vL[2])));
      pR = (gamma-1.0)*(UR(fpt,4) - (0.5*UR(fpt,0)*(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2])));
    }

    double hL = (UL(fpt,nDims+1)+pL)/UL(fpt,0);
    double hR = (UR(fpt,nDims+1)+pR)/UR(fpt,0);

    double sq_rho = sqrt(UR(fpt,0)/UL(fpt,0));
    double rrho = 1./(sq_rho+1.);

    // Roe-averaged velocity
    double usq = 0.;
    double unm = 0.;
    double vgn = 0.;
    for (int i=0; i<nDims; i++) {
      um[i] = rrho*(vL[i]+sq_rho*vR[i]);
      usq += 0.5*um[i]*um[i];
      unm += um[i]*normL(fpt,i);
    }

    // Roe-averaged enthalpy
    double hm = rrho*(hL + sq_rho*hR);

    // Roe-avereged speed of sound
    double am_sq = (gamma-1.)*(hm-usq);
    double am = sqrt(am_sq);

    // Compute Euler flux (first part)
    double rhoUnL = 0.;
    double rhoUnR = 0.;
    for (int i=0;i<nDims;i++) {
      rhoUnL += UL(fpt,i+1)*normL(fpt,i);
      rhoUnR += UR(fpt,i+1)*normL(fpt,i);
    }

    Fn(fpt,0) = rhoUnL + rhoUnR;
    Fn(fpt,1) = rhoUnL*vL[0] + rhoUnR*vR[0] + (pL+pR)*normL(fpt,0);
    Fn(fpt,2) = rhoUnL*vL[1] + rhoUnR*vR[1] + (pL+pR)*normL(fpt,1);
    Fn(fpt,3) = rhoUnL*hL   +rhoUnR*hR;

    // Compute solution difference across face
    for (int i=0;i<params->nFields;i++) {
      du[i] = UR(fpt,i)-UL(fpt,i);
    }

    // Eigenvalues
    double lambda0 = abs(unm-vgn);
    double lambdaP = abs(unm-vgn+am);
    double lambdaM = abs(unm-vgn-am);

    // Entropy fix
    double eps = 0.5*(abs(rhoUnL/UL(fpt,0)-rhoUnR/UR(fpt,0))+ abs(sqrt(gamma*pL/UL(fpt,0))-sqrt(gamma*pR/UR(fpt,0))));
    if(lambda0 < 2.*eps)
      lambda0 = 0.25*lambda0*lambda0/eps + eps;
    if(lambdaP < 2.*eps)
      lambdaP = 0.25*lambdaP*lambdaP/eps + eps;
    if(lambdaM < 2.*eps)
      lambdaM = 0.25*lambdaM*lambdaM/eps + eps;

    double a2 = 0.5*(lambdaP+lambdaM)-lambda0;
    double a3 = 0.5*(lambdaP-lambdaM)/am;
    double a1 = a2*(gamma-1.)/am_sq;
    double a4 = a3*(gamma-1.);

    double a5, a6;
    if (nDims==2) {
      a5 = usq*du[0]-um[0]*du[1]-um[1]*du[2]+du[3];
      a6 = unm*du[0]-normL(fpt,0)*du[1]-normL(fpt,1)*du[2];
    }
    else if (nDims==3) {
      a5 = usq*du[0]-um[0]*du[1]-um[1]*du[2]-um[2]*du[3]+du[4];
      a6 = unm*du[0]-normL(fpt,0)*du[1]-normL(fpt,1)*du[2]-normL(fpt,2)*du[3];
    }

    double aL1 = a1*a5 - a3*a6;
    double bL1 = a4*a5 - a2*a6;

    // Compute Euler flux (second part)
    if (nDims==2) {
      Fn(fpt,0) -= lambda0*du[0]+aL1;
      Fn(fpt,1) -= lambda0*du[1]+aL1*um[0]+bL1*normL(fpt,0);
      Fn(fpt,2) -= lambda0*du[2]+aL1*um[1]+bL1*normL(fpt,1);
      Fn(fpt,3) -= lambda0*du[3]+aL1*hm   +bL1*unm;
    }
    else if (nDims==3) {
      Fn(fpt,0) -= lambda0*du[0]+aL1;
      Fn(fpt,1) -= lambda0*du[1]+aL1*um[0]+bL1*normL(fpt,0);
      Fn(fpt,2) -= lambda0*du[2]+aL1*um[1]+bL1*normL(fpt,1);
      Fn(fpt,3) -= lambda0*du[3]+aL1*um[2]+bL1*normL(fpt,2);
      Fn(fpt,4) -= lambda0*du[4]+aL1*hm   +bL1*unm;
    }

    for (int i=0;i<params->nFields;i++) {
      Fn(fpt,i) *= 0.5;
    }
  }
}

void face::laxFriedrichsFlux(void)
{
  for (int fpt=0; fpt<nFptsL; fpt++) {
    if (params->equation == ADVECTION_DIFFUSION) {
      double uAvg = 0.5*(UL(fpt,0) + UR(fpt,0));
      double uDiff = UL(fpt,0) - UR(fpt,0);
      double vNorm = params->advectVx*normL(fpt,0) + params->advectVy*normL(fpt,1);
      if (nDims==3)
        vNorm += params->advectVz*normL(fpt,2);

      Fn(fpt,0) = vNorm*uAvg + 0.5*params->lambda*abs(vNorm)*uDiff;

      vNorm = abs(vNorm);

      if (params->motion) {
        double vgN = Vg(fpt,0) * normL(fpt,0) + Vg(fpt,1) * normL(fpt,1);
        if (nDims == 3) vgN += Vg(fpt,2) * normL(fpt,2);
        double vTot = abs(vNorm - vgN);
        vNorm = max(vTot,max(vNorm,abs(vgN)));
      }

      *waveSp[fpt] = vNorm;
    }
    else if (params->equation == NAVIER_STOKES) {
      FatalError("Lax-Friedrichs not supported for Navier-Stokes simulations - use Rusanov or Roe.");
    }
  }
}

void face::centralFluxBound(void)
{
  for (int fpt=0; fpt<nFptsL; fpt++) {
    inviscidFlux(UL[fpt],tempFL,params);
    inviscidFlux(UR[fpt],tempFR,params);

    double tempFnL[5] = {0,0,0,0,0};
    double tempFnR[5] = {0,0,0,0,0};

    // Get primitive variables
    double rhoL = UL(fpt,0);     double rhoR = UR(fpt,0);
    double uL = UL(fpt,1)/rhoL;  double uR = UR(fpt,1)/rhoR;
    double vL = UL(fpt,2)/rhoL;  double vR = UR(fpt,2)/rhoR;

    double wL, pL, vnL=0.;
    double wR, pR, vnR=0.;
    double vgn=0.;

    // Calculate pressure
    if (params->nDims==2) {
      pL = (params->gamma-1.0)*(UL(fpt,3)-rhoL*(uL*uL+vL*vL));
      pR = (params->gamma-1.0)*(UR(fpt,3)-rhoR*(uR*uR+vR*vR));
    }
    else {
      wL = UL(fpt,3)/rhoL;   wR = UR(fpt,3)/rhoR;
      pL = (params->gamma-1.0)*(UL(fpt,4)-rhoL*(uL*uL+vL*vL+wL*wL));
      pR = (params->gamma-1.0)*(UR(fpt,4)-rhoR*(uR*uR+vR*vR+wR*wR));
    }

    // Get normal fluxes, normal velocities
    for (int dim=0; dim<params->nDims; dim++) {
      vnL += normL(fpt,dim)*UL(fpt,dim+1)/rhoL;
      vnR += normL(fpt,dim)*UR(fpt,dim+1)/rhoR;
      if (params->motion)
        vgn += normL(fpt,dim)*Vg(fpt,dim);
      for (int i=0; i<params->nFields; i++) {
        tempFnL[i] += normL(fpt,dim)*tempFL[dim][i];
        tempFnR[i] += normL(fpt,dim)*tempFR[dim][i];
      }
    }

    // Get maximum eigenvalue for diffusion coefficient
    double csqL = max(params->gamma*pL/rhoL,0.0);
    double csqR = max(params->gamma*pR/rhoR,0.0);
    double eigL = std::fabs(vnL) + sqrt(csqL);
    double eigR = std::fabs(vnR) + sqrt(csqR);

    // Calculate Rusanov flux
    for (int i=0; i<params->nFields; i++) {
      Fn(fpt,i) = 0.5*(tempFnL[i]+tempFnR[i]);
    }

    // Store wave speed for calculation of allowable dt
    if (params->motion) {
      eigL = std::fabs(vnL-vgn) + sqrt(csqL);
      eigR = std::fabs(vnR-vgn) + sqrt(csqR);
    }
    *waveSp[fpt] = max(eigL,eigR);
  }
}

//! First step of the LDG flux - take a biased average of the solution
void face::ldgSolution()
{
  if (isBnd) {
    for (int fpt=0; fpt<nFptsL; fpt++)
      for(int k=0;k<nFields;k++)
        UC(fpt,k) = 0.5*(UL(fpt,k) + UR(fpt,k));
  }
  else {
    // Choosing a unique direction for the switch
    for (int fpt=0; fpt<nFptsL; fpt++) {

      double penFact = params->penFact;

      if (nDims == 2) {
        if ( normL(fpt,0)+normL(fpt,1) < 0 )
          penFact = -params->penFact;
      }
      else if (nDims == 3) {
        if (normL(fpt,0)+normL(fpt,1)+sqrt(2.)*normL(fpt,2) < 0) {
          penFact = -params->penFact;
        }
      }

      for(int k=0;k<nFields;k++)
        UC(fpt,k) = 0.5*(UL(fpt,k) + UR(fpt,k)) - penFact*(UL(fpt,k) - UR(fpt,k));
    }
  }
}

point face::getPosFpt(int fpt)
{
  return eL->getPosFpt(fptStartL+fpt);
}
