/*!
 * \file points.hpp
 * \brief Functions related to solution, flux, and quadrature points
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */
#pragma once

#include <vector>

#include "global.hpp"

//! Get the reference-domain location of the solution points for the given element & polynomial order
vector<point> getLocSpts(int eType, int order, string sptsType);

//! Get the reference-domain location of the flux points for the given element & polynomial order
vector<point> getLocFpts(int eType, int order, string sptsType);

//! Get the point locations of the requested type (i.e. Gauss, Lobatto) for the given order
vector<double> getPts1D(string ptsType, int order);

//! Get the Gauss quadrature weights for the Gauss points of the given order [2D]
vector<double> getQptWeights(int order, int nDims);

//! Get the Gauss quadrature weights for the Gauss points of the given order [1D]
vector<double> getQptWeights1D(int order);

//! Get quadrature rule (points & weights) for a tetrahedron for a given order
void getQuadRuleTet(int order, vector<point> &locQpts, vector<double> &weights);
