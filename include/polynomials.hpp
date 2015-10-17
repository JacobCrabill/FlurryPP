/*!
 * \file polynomials.hpp
 * \brief Header file for polynomial functions
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */

#include "global.hpp"

/*! Evaluate the 1D Lagrange polynomial mode based on points x_lag at point y */
double Lagrange(vector<double> &x_lag, double y, uint mode);

/*! Evaluate the first derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double dLagrange(vector<double> &x_lag, double y, uint mode);

/*! Evaluate the second derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double ddLagrange(vector<double> &x_lag, double y, uint mode);

/*! Evaluate the 1D Legendre polynomial mode number in_mode based at point location in_r*/
double Legendre(double in_r, int in_mode);

/*! Evaluate the derivative of the 1D Legendre polynomial mode number in_mode based at point location in_r*/
double dLegendre(double in_r, int in_mode);

/*! Evaluate the 2D Legendre polynomial mode number in_mode based at point location in_r*/
double Legendre2D_hierarchical(int in_mode, vector<double> in_loc, int in_basis_order);

/*! Method to calculate co-efficients of the exponential filter */
double exponential_filter(int in_mode, int inBasisOrder, double exponent);

double compute_eta(int vcjh_scheme, int order);

/*! Evaluate the divergence [dF/dxi or dG/deta] of the correction function at a point */
double dVCJH_1d(double in_r, int in_mode, int in_order, double in_eta);

// correction function bases

// Triangule Dubiner basis

//! Evaluate a mode of the [alpha,beta] Jacobi polynomial at a point
double eval_jacobi(double in_r, int in_alpha, int in_beta, int in_mode);

//! Evaluate the derivative of a normalized Jacobi polynomial
double eval_grad_jacobi(double in_r, int in_alpha, int in_beta, int in_mode);

//! Evaluate the Dubiner basis for triangles
double eval_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order);

//! Evaluate the r derivative of the Dubiner basis for triangles
double eval_dr_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order);

//! Evaluate the s derivative of the Dubiner basis for triangles
double eval_ds_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order);

//! Evaluate the gamma function for positive integers
double eval_gamma(int in_n);

//! Change triangular coordinates from natural to area
point rs_to_ab(point &rs);
