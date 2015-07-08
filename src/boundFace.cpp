/*!
 * \file boundFace.cpp
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
  bcType = myInfo.bcType;

  if (params->slipPenalty) {
    // For PID-controlled BC
    deltaU.setup(nFptsL,nDims);
    deltaUdot.setup(nFptsL,nDims);
    deltaUint.setup(nFptsL,nDims);
    deltaU.initializeToZero();
    deltaUdot.initializeToZero();
    deltaUint.initializeToZero();
    UR.initializeToZero();
  }
}

void boundFace::getRightState(void)
{
  // Set the boundary condition [store in UR]
  applyBCs();
}

void boundFace::applyBCs(void)
{
  uint nDims = params->nDims;
  for (int fpt=0; fpt<nFptsL; fpt++) {

    if (params->equation == NAVIER_STOKES) {
      // These varibles will be used to set the right state of the boundary.
      double rhoR, pR, ER, TR;
      array<double,3> vL = {0,0,0};
      array<double,3> vR = {0,0,0};
      array<double,3> vG = {0,0,0};
      array<double,3> vBound = {params->uBound,params->vBound,params->wBound};

      double gamma = params->gamma;

      /* --- Calcualte primitives on left side (interior) --- */
      double rhoL = UL(fpt,0);
      double eL = UL(fpt,nDims+1);
      for (uint i=0; i<nDims; i++)
        vL[i] = UL(fpt,i+1)/UL(fpt,0);

      double vSq = 0;
      for (uint i=0; i<nDims; i++)
        vSq += (vL[i]*vL[i]);

      // --------- For PID b.c. controller -----------
      if (UR(fpt,0)==0) {
        UR(fpt,0)= UL(fpt,0);
        UR(fpt,1) = UL(fpt,1);
        UR(fpt,2) = UL(fpt,2);
      }
      vR[0] = UR(fpt,1)/UR(fpt,0);
      vR[1] = UR(fpt,2)/UR(fpt,0);
      // ---------------------------------------------

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
        vR[2] = params->wBound;

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
      else if(bcType == SLIP_WALL || bcType == SYMMETRY) {
        // extrapolate density
        rhoR = rhoL;

        // Compute normal velocity on left side
        double vnL = 0.;
        for (uint i=0; i<nDims; i++)
          vnL += (vL[i]-vG[i])*normL(fpt,i);

        // reflect normal velocity
        if (params->slipPenalty) {
          double duOld[2];
          for (uint i=0; i<nDims; i++)
            duOld[i] = deltaU(fpt,i);

          for (uint i=0; i<nDims; i++) {
            double u_bc = vL[i] - (2.0)*vnL*normL(fpt,i);
            deltaU(fpt,i) = u_bc - vR[i];
            deltaUdot(fpt,i) = (deltaU(fpt,i) - duOld[i]) / params->dt;
            deltaUint(fpt,i)+= (deltaU(fpt,i) + duOld[i]) * params->dt/2.0;
            vR[i] += params->dt*(params->Kp*20.*deltaU(fpt,i) + params->Kd/10.*deltaUdot(fpt,i) + params->Ki*deltaUint(fpt,i));
          }
        }
        else {
          for (uint i=0; i<nDims; i++) {
            vR[i] = vL[i] - (2.0)*vnL*normL(fpt,i);
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
          vnL += vL[i]*normL(fpt,i);

        double vnBound = 0;
        for (uint i=0; i<nDims; i++)
          vnBound += vBound[i]*normL(fpt,i);

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
            vR[i] = vn_star*normL(fpt,i) + (vBound[i] - vnBound*normL(fpt,i));

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
            vR[i] = vn_star*normL(fpt,i) + (vL[i] - vnL*normL(fpt,i));

          pR = rhoR/gamma*cStar*cStar;
          vSq = 0.;
          for (uint i=0; i<nDims; i++)
            vSq += (vR[i]*vR[i]);
          ER = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
        }
      }
      else {
        cout << "Boundary Condition: " << bcType << endl;
        FatalError("Boundary condition not recognized.");
      }

      // Assign calculated values to right state
      UR(fpt,0) = rhoR;
      for (uint i=0; i<nDims; i++)
        UR(fpt,i+1) = rhoR*vR[i];
      UR(fpt,nDims+1) = ER;
    }
    else if (params->equation == ADVECTION_DIFFUSION) {

    }
  }
}


void boundFace::setRightStateFlux(void)
{
  // No right state; do nothing
}

void boundFace::setRightStateSolution(void)
{
  // No right state; do nothing
}
