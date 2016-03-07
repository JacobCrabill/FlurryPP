/*!
 * \file boundFace.cpp
 * \brief Class to handle enforcement of boundary conditions at boundary faces
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

#include "boundFace.hpp"

#include <array>

#include "flux.hpp"
#include "points.hpp"

void boundFace::setupRightState(void)
{
  // This is kinda messy, but avoids separate initialize function
  bcType = myInfo.bcType;

  if (bcType == ADIABATIC_NOSLIP) // For LDG numerical fluxes
    isBnd = 2;

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

void boundFace::getPointersRight(void)
{
  // Do nothing
}

void boundFace::getRightState(void)
{
  // Set the boundary condition [store in UR]
  applyBCs();
}

void boundFace::getRightGradient(void)
{
  // Re-Set the inviscid boundary condition [on UR]
  applyBCs();

  // Set the viscous boundary condition [on gradUR]
  applyViscousBCs();
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
      if (params->motion) {
        for (int dim=0; dim<nDims; dim++)
          vG[dim] = Vg(fpt,dim);
      }
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
      if (params->slipPenalty && UR(fpt,0)==0) {
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

      // Characteristic [Copied from HiFiLES]
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
          one_over_s = params->oneOverS;

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

void boundFace::applyViscousBCs(void)
{
  /* Apply the adiabatic-wall boundary condition to the energy gradient
   * (by removing ALL temparture gradients). TODO: remove ONLY normal comp. */
  if (bcType == ADIABATIC_NOSLIP) {
    gradUR = gradUL;

    for (int fpt=0; fpt<nFptsL; fpt++) {
      double rhovSq = 0.;
      for (int dim=0; dim<nDims; dim++)
        rhovSq += UL(fpt,dim+1)*UL(fpt,dim+1);
      double pL = (params->gamma-1.)*(UL(fpt,nDims+1) - 0.5*rhovSq/UL(fpt,0));

      // Extrapolate pressure; calculate internal energy
      double pR = pL;
      double e = pR/((params->gamma-1.)*UR(fpt,0));

      // Get velocity gradients
      matrix<double> gradVel(nDims,nDims);
      for (int dim1=0; dim1<nDims; dim1++)
        for (int dim2=0; dim2<nDims; dim2++)
          gradVel(dim1,dim2) = (gradUR[fpt](dim1,dim2+1) - gradUR[fpt](dim1,0)*UR(fpt,dim2+1)/UR(fpt,0))/UR(fpt,0);

      // Set energy gradient (set gradT = 0) (TODO: only remove dT_d[wall normal])
      double vSq = 0.;
      for (int dim=0; dim<nDims; dim++)
        vSq += UR(fpt,dim+1)*UR(fpt,dim+1)/(UR(fpt,0)*UR(fpt,0));

      for (int dim1=0; dim1<nDims; dim1++) {
        gradUR[fpt](dim1,nDims+1) = (e+0.5*vSq)*gradUR[fpt](dim1,0);
        for (int dim2=0; dim2<nDims; dim2++) {
          gradUR[fpt](dim1,nDims+1) += UR(fpt,dim2+1)*gradVel(dim2,dim1);
        }
      }
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

vector<double> boundFace::computeWallForce(void)
{
  vector<double> force = {0,0,0,0,0,0};

  if (bcType == SLIP_WALL || bcType == ADIABATIC_NOSLIP || bcType == ISOTHERMAL_NOSLIP) {
    int order;
    if (params->nDims == 2)
      order = nFptsL-1;
    else
      order = sqrt(nFptsL)-1;

    auto weights = getQptWeights1D(order);

    for (int fpt=0; fpt<nFptsL; fpt++) {
      double weight;
      if (nDims==2) {
        weight = weights[fpt];
      } else {
        int ifpt = fpt%(order+1);
        int jfpt = floor(fpt/(order+1));
        weight = weights[ifpt]*weights[jfpt];
      }

      double rho = UL(fpt,0);

      double vMagSq = 0;
      for (int dim=0; dim<nDims; dim++)
        vMagSq += UL(fpt,dim+1)*UL(fpt,dim+1)/(rho*rho);

      double p = (params->gamma-1)*(UL(fpt,nDims+1) - 0.5*rho*vMagSq);

      // Convective forces
      for (int dim=0; dim<nDims; dim++)
        force[dim] += p*normL(fpt,dim)*dAL[fpt]*weight;

      // Viscous forces
      if (params->viscous) {
        auto tau = viscousStressTensor(UL[fpt],gradUL[fpt],params);
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            force[3+dim1] -= tau(dim1,dim2)*normL(fpt,dim2)*dAL[fpt]*weight;
      }
    }
  }

  return force;
}

vector<double> boundFace::computeMassFlux(void)
{
  vector<double> flux(nFields);

  if (bcType == CHAR || bcType == SUP_IN || bcType == SUP_OUT || bcType == SUB_IN || bcType == SUB_OUT) {
    int order;
    if (params->nDims == 2)
      order = nFptsL-1;
    else
      order = sqrt(nFptsL)-1;

    auto weights = getQptWeights1D(order);

    for (int fpt=0; fpt<nFptsL; fpt++) {
      double weight;
      if (nDims==2) {
        weight = weights[fpt];
      } else {
        int ifpt = fpt%(order+1);
        int jfpt = floor(fpt/(order+1));
        weight = weights[ifpt]*weights[jfpt];
      }

      for (int k=0; k<nFields; k++)
        flux[k] += Fn(fpt,k)*weight;
    }
  }

  return flux;
}
