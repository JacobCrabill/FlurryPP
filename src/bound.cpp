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

#include "../include/bound.hpp"

#include <array>

void bound::setupBound(ele *eL, int locF_L, int bcType, int gID, input *params)
{
  int fptStartL, fptEndL;

  ID = gID;
  this->bcType = bcType;
  this->locF_L = locF_L;
  this->eL = eL;
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;

  nFptsL = eL->order+1;

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.resize(nFptsL);
  UR.setup(nFptsL,nFields);
  FL.resize(nFptsL);
  Fn.setup(nFptsL,nFields);
  normL.setup(nFptsL,nDims);
  FnL.resize(nFptsL);
  dAL.resize(nFptsL);
  detJacL.resize(nFptsL);

  // Get access to data at left element
  int fpt=0;
  for (int i=fptStartL; i<fptEndL; i++) {
    UL[fpt] = (eL->U_fpts[i]);
    FnL[fpt] = (eL->Fn_fpts[i]);
    dAL[fpt] = (eL->dA_fpts[i]);
    detJacL[fpt] = (eL->detJac_fpts[i]); // change to double**[] = &(eL->det[])

    for (int dim=0; dim<nDims; dim++)
      normL[fpt][dim] = (eL->norm_fpts[i][dim]); // change to dbl ptr

    FL[fpt].setup(nDims,nFields);
    for (int dim=0; dim<nDims; dim++)
      for (int k=0; k<nFields; k++)
        FL[fpt][dim][k] = &(eL->F_fpts[dim][i][k]);

    fpt++;
  }

  // Setup a temporary flux-storage vector for later use
  tempFL.setup(nDims,nFields);
  tempFR.setup(nDims,nFields);
  tempUL.resize(nFields);
}

void bound::calcInviscidFlux()
{
  for (int i=0; i<nFptsL; i++) {
    // Set the boundary condition [store in UC]
    applyBCs(UL[i],UR[i],normL[i]);

    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      centralFlux(UL[i], UR[i], normL[i], Fn[i], params);
    }
    else if (params->equation == NAVIER_STOKES) {
      if (params->riemann_type==0) {
        inviscidFlux(UL[i],tempFL,params);
        inviscidFlux(UR[i],tempFR,params);
        centralFlux(tempFL, tempFR, normL[i], Fn[i], params);
      }
      else if (params->riemann_type==1) {
        roeFlux(UL[i], UR[i], normL[i], Fn[i], params);
      }
    }

    // Calculate difference between discontinuous & common normal flux, and store in ele
    // (Each ele needs only the difference, not the actual common value, for the correction)
    // Need dAL/R to transform normal flux back to reference space
    for (int j=0; j<nFields; j++) {
      FnL[i][j] = Fn[i][j]*dAL[i];
    }
  }
}

void bound::calcViscousFlux()
{

}

void bound::applyBCs(const double* uL, double* uR, const double *norm)
{
  uint nDims = params->nDims;

  if (params->equation == NAVIER_STOKES) {
    // These varibles will be used to set the right state of the boundary.
    double rhoR, pR, eR, TR;
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

    // --------- TESTING -----------
    if (uR[0]==0) {
      uR[0] = uL[0];
      uR[1] = uL[1];
      uR[2] = uL[2];
    }
    vR[0] = uR[1]/uR[0];
    vR[1] = uR[2]/uR[0];
    // --------- TESTING -----------

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
      eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
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

      eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }

    // Subsonic inflow characteristic
    // there is one outgoing characteristic (u-c), therefore we can specify
    // all but one state variable at the inlet. The outgoing Riemann invariant
    // provides the final piece of info. Adapted from an implementation in
    // SU2.
    else if(bcType == SUB_IN_CHAR) {
      /*
        double vR;
        double cL, cR_sq, c_total_sq;
        double R_plus, h_total;
        double aa, bb, cc, dd;
        double Mach_sq, alpha;

        // Specify Inlet conditions
        double p_total_bound = bdy_params[9];
        double T_total_bound = bdy_params[10];
        double *n_free_stream = &bdy_params[11];

        // Compute normal velocity on left side
        vnL = 0.;
        for (int i=0; i<nDims; i++)
          vnL += vL[i]*norm[i];

        // Compute speed of sound
        cL = sqrt(gamma*pL/rhoL);

        // Extrapolate Riemann invariant
        R_plus = vnL + 2.0*cL/(gamma-1.0);

        // Specify total enthalpy
        h_total = gamma*R_ref/(gamma-1.0)*T_total_bound;

        // Compute total speed of sound squared
        vSq = 0.;
        for (int i=0; i<nDims; i++)
          vSq += vL[i]*vL[i];

        c_total_sq = (gamma-1.0)*(h_total - (eL/rhoL + pL/rhoL) + 0.5*vSq) + cL*cL;

        // Dot product of normal flow velocity
        alpha = 0.;
        for (int i=0; i<nDims; i++)
          alpha += norm[i]*n_free_stream[i];

        // Coefficients of quadratic equation
        aa = 1.0 + 0.5*(gamma-1.0)*alpha*alpha;
        bb = -(gamma-1.0)*alpha*R_plus;
        cc = 0.5*(gamma-1.0)*R_plus*R_plus - 2.0*c_total_sq/(gamma-1.0);

        // Solve quadratic equation for velocity on right side
        // (Note: largest value will always be the positive root)
        // (Note: Will be set to zero if NaN)
        dd = bb*bb - 4.0*aa*cc;
        dd = sqrt(max(dd, 0.0));
        vR = (-bb + dd)/(2.0*aa);
        vR = max(vR, 0.0);
        vSq = vR*vR;

        // Compute speed of sound
        cR_sq = c_total_sq - 0.5*(gamma-1.0)*vSq;

        // Compute Mach number (cutoff at Mach = 1.0)
        Mach_sq = vSq/(cR_sq);
        Mach_sq = min(Mach_sq, 1.0);
        vSq = Mach_sq*cR_sq;
        vR = sqrt(vSq);
        cR_sq = c_total_sq - 0.5*(gamma-1.0)*vSq;

        // Compute velocity (based on free stream direction)
        for (int i=0; i<nDims; i++)
          vR[i] = vR*n_free_stream[i];

        // Compute temperature
        TR = cR_sq/(gamma*R_ref);

        // Compute pressure
        pR = p_total_bound*pow(TR/T_total_bound, gamma/(gamma-1.0));

        // Compute density
        rhoR = pR/(R_ref*TR);

        // Compute energy
        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
        */
    }

    // Subsonic outflow characteristic
    // there is one incoming characteristic, therefore one variable can be
    // specified (back pressure) and is used to update the conservative
    // variables. Compute the entropy and the acoustic Riemann variable.
    // These invariants, as well as the tangential velocity components,
    // are extrapolated. Adapted from an implementation in SU2.
    else if(bcType == SUB_OUT_CHAR) {
      /*
        double cL, cR;
        double R_plus, s;
        double vnR;

        // Compute normal velocity on left side
        vnL = 0.;
        for (int i=0; i<nDims; i++)
          vnL += vL[i]*norm[i];

        // Compute speed of sound
        cL = sqrt(gamma*pL/rhoL);

        // Extrapolate Riemann invariant
        R_plus = vnL + 2.0*cL/(gamma-1.0);

        // Extrapolate entropy
        s = pL/pow(rhoL,gamma);

        // fix pressure on the right side
        pR = pBound;

        // Compute density
        rhoR = pow(pR/s, 1.0/gamma);

        // Compute speed of sound
        cR = sqrt(gamma*pR/rhoR);

        // Compute normal velocity
        vnR = R_plus - 2.0*cR/(gamma-1.0);

        // Compute velocity and energy
        vSq = 0.;
        for (int i=0; i<nDims; i++) {
          vR[i] = vL[i] + (vnR - vnL)*norm[i];
          vSq += (vR[i]*vR[i]);
        }
        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
        */
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

      eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
    }


    // Supersonic outflow
    else if(bcType == SUP_OUT) {
      // extrapolate density, velocity, energy
      rhoR = rhoL;
      for (uint i=0; i<nDims; i++)
        vR[i] = vL[i];
      eR = eL;
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
        for (uint i=0; i<nDims; i++) {
          vR[i] = vR[i] - params->beta*params->dt*(vR[i] - (vL[i] - (2.0)*vnL*norm[i]));
          //vR[i] = vL[i] - (2.0-(double)5000./(params->iter+2500.))*vnL*norm[i];
        }
      }
      else {
        for (uint i=0; i<nDims; i++) {
          vR[i] = vL[i] - (2.0)*vnL*norm[i];
        }
      }

      // extrapolate energy
      eR = eL;
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

      eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
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

      eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
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
        eR = rhoR*h_free_stream - pR;
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
        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
      }
    }

    // Assign calculated values to right state
    uR[0] = rhoR;
    for (uint i=0; i<nDims; i++)
      uR[i+1] = rhoR*vR[i];
    uR[nDims+1] = eR;
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    // Trivial Dirichlet
    /*if(bdy_type==50)
      {
        uR[0]=0.0;
      }*/
  }
}
