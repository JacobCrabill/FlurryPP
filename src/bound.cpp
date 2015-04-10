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

void bound::setupBound(ele *eL, int locF_L, int bcType, int gID)
{
  int fptStartL, fptEndL;

  ID = gID;
  this->bcType = bcType;
  this->locF_L = locF_L;

  nDims = params->nDims;
  nFields = params->nFields;

  nFptsL = eL->order+1;

  /* --- For 1D faces [line segments] only - find first/last ID of fpts; reverse
   * the order on the 'right' face so they match up --- */
  fptStartL = (locF_L*(nFptsL));
  fptEndL = (locF_L*(nFptsL)) + nFptsL;

  UL.resize(nFptsL);
  FL.resize(nFptsL);
  normL.setup(nFptsL,nDims);
  detJacL.resize(nFptsL);

  // Get access to data at left element
  int fpt=0;
  for (int i=fptStartL; i<fptEndL; i++) {
    UL[fpt] = (eL->U_fpts[i]);
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
  tempUL.resize(nFields);
}

void bound::calcInviscidFlux()
{
  for (int i=0; i<nFptsL; i++) {
    // Set the boundary condition [store in UC]
    applyBCs(UL[i],&UC[i][0],normL[i]);

    // Calcualte discontinuous inviscid flux at flux points
    inviscidFlux(UL[i],tempFL, params);

    // Calculate common inviscid flux at flux points
    if (params->equation == ADVECTION_DIFFUSION) {
      upwindFlux(UL[i], &UC[i][0], normL[i], Fn[i], params);
    }
//    else if (params->equation == NAVIER_STOKES) {
//      if (params->riemann_type==0) {
//        rusanovFlux(*UL[i], *UC[i], *FL[i], *FR[i], normL[i], *Fn[i], params);
//      }
//      else if (params->riemann_type==1) {
//        roeFlux(*UL[i], *UR[i], normL[i], *Fn[i], params);
//      }
//    }
  }
}

void bound::calcViscousFlux()
{

}

void bound::applyBCs(double* uL, double* uR, double* norm)
{
  uint nDims = params->nDims;

  for (int fpt=0; fpt<nFptsL; fpt++) {
    if (params->equation == NAVIER_STOKES) {
      double rhoL, rhoR, pL, pR, eL, eR, vSq;
      array<double,3> vL, vR;

      double gamma = params->gamma;

      /* --- Calcualte primitives on left side (interior) --- */
      rhoL = uL[0];
      eL = uL[nDims+1];
      for (int i=0; i<nDims; i++)
        vL[i] = uL[i+1]/uL[0];

      vSq = 0;
      for (int i=0; i<nDims; i++)
        vSq += (vL[i]*vL[i]);

      pL = (gamma-1.0)*(eL - 0.5*rhoL*vSq);

      // Subsonic inflow simple (free pressure) //CONSIDER DELETING
/*      if(bcType == 1) {
        // fix density and velocity
        rhoR = rho_bound;
        for (int i=0; i<nDims; i++)
          vR[i] = v_bound[i];

        // extrapolate pressure
        pR = pL;

        // compute energy
        vSq = 0;
        for (int i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);
        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
      }

      // Subsonic outflow simple (fixed pressure) //CONSIDER DELETING
      else if(bcType == 2) {
        // extrapolate density and velocity
        rhoR = rhoL;
        for (int i=0; i<nDims; i++)
          vR[i] = vL[i];

        // fix pressure
        pR = p_bound;

        // compute energy
        vSq = 0.;
        for (int i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);

        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
      }

      // Subsonic inflow characteristic
      // there is one outgoing characteristic (u-c), therefore we can specify
      // all but one state variable at the inlet. The outgoing Riemann invariant
      // provides the final piece of info. Adapted from an implementation in
      // SU2.
      else if(bdy_type == 3) {
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

        // SA model
        if (run_input.turb_model == 1) {
          // set turbulent eddy viscosity
          double mu_tilde_inf = bdy_params[14];
          uR[nDims+2] = mu_tilde_inf;
        }
      }

      // Subsonic outflow characteristic
      // there is one incoming characteristic, therefore one variable can be
      // specified (back pressure) and is used to update the conservative
      // variables. Compute the entropy and the acoustic Riemann variable.
      // These invariants, as well as the tangential velocity components,
      // are extrapolated. Adapted from an implementation in SU2.
      else if(bdy_type == 4) {
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
        pR = p_bound;

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

      }

      // Supersonic inflow
      else if(bdy_type == 5) {
        // fix density and velocity
        rhoR = rho_bound;
        for (int i=0; i<nDims; i++)
          vR[i] = v_bound[i];

        // fix pressure
        pR = p_bound;

        // compute energy
        vSq = 0.;
        for (int i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);
        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
      }


      // Supersonic outflow
      else if(bdy_type == 6) {
        // extrapolate density, velocity, energy
        rhoR = rhoL;
        for (int i=0; i<nDims; i++)
          vR[i] = vL[i];
        eR = eL;
      }

      // Slip wall
      else if(bdy_type == 7) {
        // extrapolate density
        rhoR = rhoL;

        // Compute normal velocity on left side
        vnL = 0.;
        for (int i=0; i<nDims; i++)
          vnL += (vL[i]-v_g[i])*norm[i];

        // reflect normal velocity
        for (int i=0; i<nDims; i++)
          vR[i] = vL[i] - 2.0*vnL*norm[i];

        // extrapolate energy
        eR = eL;
      }

      // Isothermal, no-slip wall (fixed)
      else if(bdy_type == 11) {
        // Set state for the right side
        // extrapolate pressure
        pR = pL;

        // isothermal temperature
        TR = T_wall;

        // density
        rhoR = pR/(R_ref*TR);

        // no-slip
        for (int i=0; i<nDims; i++)
          vR[i] = v_g[i];

        // energy
        vSq = 0.;
        for (int i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);

        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;

        // SA model
        if (run_input.turb_model == 1) {
          // zero turbulent eddy viscosity at the wall
          uR[nDims+2] = 0.0;
        }
      }

      // Adiabatic, no-slip wall (fixed)
      else if(bdy_type == 12) {
        // extrapolate density
        rhoR = rhoL; // only useful part

        // extrapolate pressure
        pR = pL;

        // no-slip
        for (int i=0; i<nDims; i++)
          vR[i] = v_g[i];

        // energy
        vSq = 0.;
        for (int i=0; i<nDims; i++)
          vSq += (vR[i]*vR[i]);

        eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;

        // SA model
        if (run_input.turb_model == 1) {
          // zero turbulent eddy viscosity at the wall
          uR[nDims+2] = 0.0;
        }
      }

      // Characteristic
      else if (bdy_type == 15) {
        double c_star;
        double vn_star;
        double vn_bound;
        double r_plus,r_minus;

        double one_over_s;
        double h_free_stream;

        // Compute normal velocity on left side
        vnL = 0.;
        for (int i=0; i<nDims; i++)
          vnL += vL[i]*norm[i];

        vn_bound = 0;
        for (int i=0; i<nDims; i++)
          vn_bound += v_bound[i]*norm[i];

        r_plus  = vnL + 2./(gamma-1.)*sqrt(gamma*pL/rhoL);
        r_minus = vn_bound - 2./(gamma-1.)*sqrt(gamma*p_bound/rho_bound);

        c_star = 0.25*(gamma-1.)*(r_plus-r_minus);
        vn_star = 0.5*(r_plus+r_minus);

        // Inflow
        if (vnL<0) {
          // HACK
          one_over_s = pow(rho_bound,gamma)/p_bound;

          // freestream total enthalpy
          vSq = 0.;
          for (int i=0;i<nDims;i++)
            vSq += v_bound[i]*v_bound[i];
          h_free_stream = gamma/(gamma-1.)*p_bound/rho_bound + 0.5*vSq;

          rhoR = pow(1./gamma*(one_over_s*c_star*c_star),1./(gamma-1.));

          // Compute velocity on the right side
          for (int i=0; i<nDims; i++)
            vR[i] = vn_star*norm[i] + (v_bound[i] - vn_bound*norm[i]);

          pR = rhoR/gamma*c_star*c_star;
          eR = rhoR*h_free_stream - pR;
        }
        // Outflow
        else {
          one_over_s = pow(rhoL,gamma)/pL;

          // freestream total enthalpy
          rhoR = pow(1./gamma*(one_over_s*c_star*c_star), 1./(gamma-1.));

          // Compute velocity on the right side
          for (int i=0; i<nDims; i++)
            vR[i] = vn_star*norm[i] + (vL[i] - vnL*norm[i]);

          pR = rhoR/gamma*c_star*c_star;
          vSq = 0.;
          for (int i=0; i<nDims; i++)
            vSq += (vR[i]*vR[i]);
          eR = (pR/(gamma-1.0)) + 0.5*rhoR*vSq;
        }
      }
*/
    }
    else if (params->equation == ADVECTION_DIFFUSION) {
      // Trivial Dirichlet
      /*if(bdy_type==50)
      {
        uR[0]=0.0;
      }*/
    }
  }
}
