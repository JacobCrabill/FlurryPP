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
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/face.hpp"

#include "../include/flux.hpp"
#include "../include/ele.hpp"

void face::initialize(ele *eL, ele *eR, int locF_L, const vector<int> &rightParams, int gID, input *params)
{
  ID = gID;

  this->locF_L = locF_L;
  this->eL = eL;
  this->eR = eR;
  this->params = params;

  // Note: this is locF_R for internal faces, bcType for boundary faces, and right ID for mpi faces
  this->rightParams = rightParams;

  nDims = params->nDims;
  nFields = params->nFields;

  // Setup temporary vectors for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);

  // Needed for MPI faces to separate communication from flux calculation
  isMPI = 0;
}

void face::setupFace(void)
{
  nFptsL = eL->order+1;

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.setup(nFptsL,nFields);
  UR.setup(nFptsL,nFields);
  FnL.resize(nFptsL);
  Fn.setup(nFptsL,nFields);
  normL.setup(nFptsL,nDims);
  dAL.resize(nFptsL);
  detJacL.resize(nFptsL);
  waveSp.resize(nFptsL);

  if (params->viscous) {
    UC.setup(nFptsL,nFields);
    UcL.resize(nFptsL);
    // just a placeholder.  Need to properly size/reorder dimensions later.
    gradUL.resize(nFptsL);
    gradUR.resize(nFptsL);
    for (auto &dU:gradUL) dU.setup(nDims,nFields);
    for (auto &dU:gradUR) dU.setup(nDims,nFields);
  }

  // Get access to data at left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    FnL[fpt] = (eL->Fn_fpts[i]);
    waveSp[fpt] = &(eL->waveSp_fpts[i]);

    if (params->viscous) {
      UcL[fpt] = (eL->Uc_fpts[i]);
    }

    fpt++;
  }

  this->setupRightState();
}

void face::getLeftState()
{
  // Get data from left element
  int fpt = 0;
  for (int i=fptStartL; i<fptEndL; i++) {
    for (int j=0; j<nFields; j++) {
      UL(fpt,j) = (eL->U_fpts(i,j));
    }

    // For dynamic grids, need to update geometry-related data
    if ((params->iter == params->initIter+1) || (params->motion != 0)) {
      for (int dim=0; dim<nDims; dim++) {
        normL(fpt,dim) = (eL->norm_fpts(i,dim));
      }
      dAL[fpt] = (eL->dA_fpts[i]);
      detJacL[fpt] = (eL->detJac_fpts[i]);
    }

    if (params->viscous) {
      for (int j=0; j<nFields; j++)
        for (int dim=0; dim<nDims; dim++)
          gradUL[fpt](dim,j) = (eL->dU_fpts[dim](i,j));
    }

    fpt++;
  }
}

void face::calcInviscidFlux(void)
{
  if (!isMPI)
    getLeftState();  // Idea: instead of using ptrs to ele data, do copy?
  this->getRightState(); // <-- makes this more general for all face types, and allows face memory to be contiguous

  // Calculate common inviscid flux at flux points
  if (params->equation == ADVECTION_DIFFUSION) {
    laxFriedrichsFlux();
  }
  else if (params->equation == NAVIER_STOKES) {
    if (params->riemann_type==0) {
      rusanovFlux();
    }
    else if (params->riemann_type==1) {
      roeFlux();
    }
  }

  // Transform normal flux using edge Jacobian and put into ele's memory
  for (int i=0; i<nFptsL; i++) {
    for (int j=0; j<nFields; j++)
      FnL[i][j] =  Fn(i,j)*dAL[i];

    *waveSp[i] /= dAL[i];
  }

  if (params->viscous) {

    ldgSolution();

    // Still assuming nFptsL = nFptsR
    for (int i=0; i<nFptsL; i++) {
      for (int j=0; j<nFields; j++) {
        UC(i,j) = 0.5*( UL(i,j) + UR(i,j) );
        UcL[i][j] = UC(i,j);
      }
    }
  }

  this->setRightStateFlux();

  if (params->viscous)
    this->setRightStateSolution();
}

void face::calcViscousFlux(void)
{
  if (params->equation == NAVIER_STOKES) {
    for (int fpt=0; fpt<nFptsL; fpt++) {
      // Calculate discontinuous viscous flux at flux points
      viscousFlux(UL[fpt], gradUL[fpt], tempFL, params);
      viscousFlux(UR[fpt], gradUR[fpt], tempFR, params);

      // Calculte common viscous flux at flux points
      //ldgFlux(UL[i], UR[i], tempFL, tempFL, Fn[i], params);

      double penFact = params->penFact;
      if ( normL(fpt,0)+normL(fpt,1) < 0 )
        penFact = -params->penFact;

      double normX = normL(fpt,0);
      double normY = normL(fpt,1);

      matrix<double> Fc(nDims,nFields);
      for(int k=0; k<nFields; k++) {
        Fc(0,k) = 0.5*(tempFL(0,k) + tempFR(0,k)) + penFact*normX*( normX*(tempFL(0,k) - tempFR(0,k)) + normY*(tempFL(1,k) - tempFR(1,k)) ) + params->tau*normX*(UL(fpt,k) - UR(fpt,k));
        Fc(1,k) = 0.5*(tempFL(1,k) + tempFR(1,k)) + penFact*normY*( normX*(tempFL(0,k) - tempFR(0,k)) + normY*(tempFL(1,k) - tempFR(1,k)) ) + params->tau*normY*(UL(fpt,k) - UR(fpt,k));
      }

      // calculate normal flux from discontinuous solution at flux points
      for(int k=0; k<nFields; k++)
        Fn(fpt,k) += Fc(0,k)*normL(fpt,0) + Fc(1,k)*normL(fpt,1);

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

    // Calculate pressure
    if (params->nDims==2) {
      pL = (params->gamma-1.0)*(UL(fpt,3)-rhoL*(uL*uL+vL*vL));
      pR = (params->gamma-1.0)*(UR(fpt,3)-rhoR*(uR*uR+vR*vR));
    }
    else {
      wL = UL(fpt,3)/rhoL;   wR = UR(fpt,3)/rhoR;
      pL = (params->gamma-1.0)*(UL(fpt,3)-rhoL*(uL*uL+vL*vL+wL*wL));
      pR = (params->gamma-1.0)*(UR(fpt,3)-rhoR*(uR*uR+vR*vR+wR*wR));
    }

    // Get normal fluxes, normal velocities
    for (int j=0; j<params->nDims; j++) {
      vnL += normL(fpt,j)*UL(fpt,j+1)/rhoL;
      vnR += normL(fpt,j)*UR(fpt,j+1)/rhoR;
      for (int i=0; i<params->nFields; i++) {
        tempFnL[i] += normL(fpt,j)*tempFL(j,i);
        tempFnR[i] += normL(fpt,j)*tempFR(j,i);
      }
    }

    // Get maximum eigenvalue for diffusion coefficient
    double csqL = max(params->gamma*pL/rhoL,0.0);
    double csqR = max(params->gamma*pR/rhoR,0.0);
    double eigL = std::fabs(vnL) + sqrt(csqL);
    double eigR = std::fabs(vnR) + sqrt(csqR);
    *waveSp[fpt] = max(eigL,eigR);

    // Calculate Rusanov flux
    for (int i=0; i<params->nFields; i++) {
      Fn(fpt,i) = 0.5*(tempFnL[i]+tempFnR[i] - (*waveSp[fpt])*(UR(fpt,i)-UL(fpt,i)));
    }
  }
}

void face::roeFlux(void)
{
  double gamma = params->gamma;
  int nDims = params->nDims;

  if (nDims == 3) FatalError("Roe not implemented in 3D");

  for (int fpt=0; fpt<nFptsL; fpt++) {
    array<double,3> vL, vR, um;
    array<double,5> du;

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
    }
    else if (params->equation == NAVIER_STOKES) {
      FatalError("Lax-Friedrichs not supported for Navier-Stokes simulations - use Rusanov or Roe.");
    }
  }
}

//! First step of the LDG flux - take a biased average of the solution
void face::ldgSolution()
{
  // Choosing a unique direction for the switch
  for (int fpt=0; fpt<nFptsL; fpt++) {

    double penFact = params->penFact;
    if ( normL(fpt,0)+normL(fpt,1) < 0 )
      penFact = -params->penFact;

    for(int k=0;k<nFields;k++)
      UC(fpt,k) = 0.5*(UL(fpt,k) + UR(fpt,k)) - penFact*(UL(fpt,k) - UR(fpt,k));

  }
}
