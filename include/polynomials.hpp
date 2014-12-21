/*!
 * \file polynomials.hpp
 * \brief Header file for polynomial functions
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

#include "global.hpp"

/*! Evaluate the 1D Lagrange polynomial mode based on points x_lag at point y */
double Lagrange(vector<double> &x_lag, double &y, int &mode);

/*! Evaluate the first derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double dLagrange(vector<double> &x_lag, double &y, int &mode);

/*! Evaluate the second derivative of the 1D Lagrange polynomial mode based on points x_lag at point y */
double ddLagrange(vector<double> &x_lag, double &y, int &mode);

// correction function bases

// Triangule Dubiner basis

// Bilinear / Biquadratic / etc. shape functions

//! Shape function for linear quad (TODO: Generalize to N-noded quad)
void shape_quad(point &in_rs, vector<double> &out_shape);

//! Derivative of shape functions for linear quad
void dshape_quad(point &in_rs, vector<vector<double>> &out_dshape);

//! Shape function for linear triangle (TODO: Generalize to N-noded tri)
void shape_tri(point &in_rs, vector<double> &out_shape);

//! Derivative of shape functions for linear triangle
void dshape_tri(point &in_rs, vector<vector<double>> &out_dshape);

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
