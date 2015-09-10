/*!
 * \file flux.hpp
 * \brief Header file flux-calculation functions
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
#pragma once

#include "global.hpp"

#include "input.hpp"
#include "matrix.hpp"

/*! Calculate the inviscid portion of the Euler or Navier-Stokes flux vector at a point */
void inviscidFlux(double* U, matrix<double> &F, input *params);

/*! Calculate the viscous portion of the Navier-Stokes flux vector at a point */
void viscousFlux(double *U, matrix<double> &gradU, matrix<double> &Fvis, input *params);

/*! Calculate the viscous shear-stress tensor at a point */
matrix<double> viscousStressTensor(double* U, matrix<double> &gradU, input *params);

/*! Calculate the viscous flux for Advection-Diffusion */
void viscousFluxAD(matrix<double> &gradU, matrix<double> &Fvis, input *params);

/*! Calculate the common inviscid flux at a point using Roe's method */
void roeFlux(double* uL, double* uR, double *norm, double *Fn, input *params);

/*! Calculate the common inviscid flux at a point using the Rusanov scalar-diffusion method */
void rusanovFlux(double* UL, double* UR, matrix<double> &FL, matrix<double> &FR, double *norm, double *Fn, double *waveSp, input *params);

/*! Simple central-difference flux (For advection problems) */
void centralFlux(double* uL, double* uR, double *norm, double *Fn, input *params);

/*! Simple central-difference flux (For Navier-Stokes problems) */
void centralFlux(matrix<double> &FL, matrix<double> &FR, double* norm, double* Fn, input *params);

/*! Simple upwinded flux (primarily for advection problems) */
void upwindFlux(double* uL, double* uR, double *norm, double *Fn, input *params);

/*! Lax-Friedrichs flux (advection-diffusion) */
void laxFriedrichsFlux(double* uL, double* uR, double *norm, double *Fn, input *params);

/*! Calculate the common viscous flux at a point using the LDG penalty method */
void ldgFlux(double* uL, double* uR, matrix<double> &gradU_L, matrix<double> &gradU_R, double *Fn, input *params);
