/*!
 * \file bound.cpp
 * \brief Class to handle enforcement of boundary conditions at boundary faces
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

#include "../include/boundFace.hpp"

#include <array>

void boundFace::setupRightState(void)
{
  // This is kinda messy, but avoids separate initialize function
  bcType = rightParam;

  deltaU.setup(nFptsL,nDims);
  deltaUdot.setup(nFptsL,nDims);
  deltaUint.setup(nFptsL,nDims);
  deltaU.initializeToZero();
  deltaUdot.initializeToZero();
  deltaUint.initializeToZero();
  UR.initializeToZero();
}

void boundFace::getRightState(void)
{
  // Set the boundary condition [store in UR]
  for (int i=0; i<nFptsL; i++) {
    applyBCs(UL[i],UR[i],normL[i],i);
  }
}

void boundFace::applyBCs(const double* uL, double* uR, const double *norm, int fpt)
{
  uint nDims = params->nDims;

  if (params->equation == NAVIER_STOKES) {
    // These varibles will be used to set the right state of the boundary.
    double rhoR, pR, ER, TR;
    array<double,3> vL = {0,0,0};
    array<double,3> vR = {0,0,0};
    array<double,3> vG = {0,0,0};
    array<double,3> vBound = {params->uBound,params->vBound,params->wBound};

    double gamma = params->gamma;

    /* --- Calcualte primitives on left side (interior) --- */
    double rhoL = uL[0];
    double eL = uL[nDims+1];
    for (uint i=0; i<nDims; i++)
      vL[i] = uL[i+1]/uL[0];

    double vSq = 0;
    for (uint i=0; i<nDims; i++)
      vSq += (vL[i]*vL[i]);

    // --------- AA222 -----------
    if (uR[0]==0) {
      uR[0] = uL[0];
      uR[1] = uL[1];
      uR[2] = uL[2];
    }
    vR[0] = uR[1]/uR[0];
    vR[1] = uR[2]/uR[0];
    // --------- AA222 -----------

    double pL = (gamma-1.0)*(eL - 0.5*rhoL*vSq);

    // Subsonic inflow simple (free pressure) //CONSIDER DELETING
    if(bcType == SUB_IN) {
      // fix density and velocity
      rhoR = params->rhoBound;
      for (uint i=0; i<nDims; i++)
        vR[i] = vBound[i];

      // extrapolate pressure
      pR = pL;

      // compute energy
      vSq = 0;
      for (uint i=0; i<nDims; i++)
        vSq += (vR[i]*vR[i]);
      ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }

    // Subsonic outflow simple (fixed pressure) //CONSIDER DELETING
    else if(bcType == SUB_OUT) {
      // extrapolate density and velocity
      rhoR = rhoL;
      for (uint i=0; i<nDims; i++)
        vR[i] = vL[i];

      // fix pressure
      pR = params->pBound;

      // compute energy
      vSq = 0.;
      for (uint i=0; i<nDims; i++)
        vSq += (vR[i]*vR[i]);

      ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }

    // Supersonic inflow
    else if(bcType == SUP_IN) {
      // fix density and velocity
      rhoR = params->rhoBound;
      vR[0] = params->uBound;
      vR[1] = params->vBound;

      // fix pressure
      pR = params->pBound;

      // compute energy
      vSq = 0.;
      for (uint i=0; i<nDims; i++)
        vSq += (vR[i]*vR[i]);

      ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }


    // Supersonic outflow
    else if(bcType == SUP_OUT) {
      // extrapolate density, velocity, energy
      rhoR = rhoL;
      for (uint i=0; i<nDims; i++)
        vR[i] = vL[i];
      ER = eL;
    }

    // Slip wall
    else if(bcType == SLIP_WALL) {
      // extrapolate density
      rhoR = rhoL;

      // Compute normal velocity on left side
      double vnL = 0.;
      for (uint i=0; i<nDims; i++)
        vnL += (vL[i]-vG[i])*norm[i];

      // reflect normal velocity
      if (params->slipPenalty) {
        double duOld[2];
        for (uint i=0; i<nDims; i++)
          duOld[i] = deltaU(fpt,i);

        for (uint i=0; i<nDims; i++) {
          double u_bc = vL[i] - (2.0)*vnL*norm[i];
          deltaU(fpt,i) = u_bc - vR[i];
          deltaUdot(fpt,i) = (deltaU(fpt,i) - duOld[i]) / params->dt;
          deltaUint(fpt,i)+= (deltaU(fpt,i) + duOld[i]) * params->dt/2.0;
          vR[i] += params->dt*(params->Kp*20.*deltaU(fpt,i) + params->Kd/10.*deltaUdot(fpt,i) + params->Ki*deltaUint(fpt,i));
        }
      }
      else {
        for (uint i=0; i<nDims; i++) {
          vR[i] = vL[i] - (2.0)*vnL*norm[i];
        }
      }

      // extrapolate energy
      ER = eL;
    }

    // Isothermal, no-slip wall (fixed)
    else if(bcType == ISOTHERMAL_NOSLIP) {
      // Set state for the right side
      // extrapolate pressure
      pR = pL;

      // isothermal temperature
      TR = params->TWall;

      // density
      rhoR = pR/(params->RGas*TR);

      // no-slip
      for (uint i=0; i<nDims; i++)
        vR[i] = vG[i];

      // energy
      vSq = 0.;
      for (uint i=0; i<nDims; i++)
        vSq += (vR[i]*vR[i]);

      ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }

    // Adiabatic, no-slip wall (fixed)
    else if(bcType == ADIABATIC_NOSLIP) {
      // extrapolate density
      rhoR = rhoL; // only useful part

      // extrapolate pressure
      pR = pL;

      // no-slip
      for (uint i=0; i<nDims; i++)
        vR[i] = vG[i];

      // energy
      vSq = 0.;
      for (uint i=0; i<nDims; i++)
        vSq += (vR[i]*vR[i]);

      ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }

    // Characteristic
    else if (bcType == CHAR) {
      double one_over_s;
      double h_free_stream;

      // Compute normal velocity on left side
      double vnL = 0.;
      for (uint i=0; i<nDims; i++)
        vnL += vL[i]*norm[i];

      double vnBound = 0;
      for (uint i=0; i<nDims; i++)
        vnBound += vBound[i]*norm[i];

      double r_plus  = vnL + 2./(gamma-1.)*sqrt(gamma*pL/rhoL);
      double r_minus = vnBound - 2./(gamma-1.)*sqrt(gamma*params->pBound/params->rhoBound);

      double cStar = 0.25*(gamma-1.)*(r_plus-r_minus);
      double vn_star = 0.5*(r_plus+r_minus);

      // Inflow
      if (vnL<0) {
        // HACK
        one_over_s = pow(params->rhoBound,gamma)/params->pBound;

        // freestream total enthalpy
        vSq = 0.;
        for (uint i=0;i<nDims;i++)
          vSq += vBound[i]*vBound[i];
        h_free_stream = gamma/(gamma-1.)*params->pBound/params->rhoBound + 0.5*vSq;

        rhoR = pow(1./gamma*(one_over_s*cStar*cStar),1./(gamma-1.));

        // Compute velocity on the right side
        for (uint i=0; i<nDims; i++)
          vR[i] = vn_star*norm[i] + (vBound[i] - vnBound*norm[i]);

        pR = rhoR/gamma*cStar*cStar;
        ER = rhoR*h_free_stream - pR;
      }
      // Outflow
      else {
        one_over_s = pow(rhoL,gamma)/pL;

        // freestream total enthalpy
        rhoR = pow(1./gamma*(one_over_s*cStar*cStar), 1./(gamma-1.));

        // Compute velocity on the right side
        for (uint i=0; i<nDims; i++)
          vR[i] = vn_star*norm[i] + (vL[i] - vnL*norm[i]);

        pR = rhoR/gamma*cStar*cStar;
        vSq = 0.;
        for (uint i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);
        ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
      }
    }
    else {
      FatalError("Boundary condition not recognized.");
    }

    // Assign calculated values to right state
    uR[0] = rhoR;
    for (uint i=0; i<nDims; i++)
      uR[i+1] = rhoR*vR[i];
    uR[nDims+1] = ER;
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    // Trivial Dirichlet
    /*if(bdy_type==50)
      {
        uR[0]=0.0;
      }*/
  }

  maxDU = std::max(maxDU,std::abs(deltaU(fpt,0))+std::abs(deltaU(fpt,1)));

}


void boundFace::setRightState(void)
{
  // No right state; do nothing
}
