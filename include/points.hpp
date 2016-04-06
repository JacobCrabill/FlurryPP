/*!
 * \file points.hpp
 * \brief Functions related to solution, flux, and quadrature points
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
#pragma once

#include <vector>

#include "global.hpp"

//! Get the reference-domain location of the solution points for the given element & polynomial order
vector<point> getLocSpts(int eType, int order, string sptsType);

//! Get the reference-domain location of the flux points for the given element & polynomial order
vector<point> getLocFpts(int eType, int order, string sptsType);

//! Get the reference-domain location of the plot points for the given element & polynomial order
vector<point> getLocPpts(int eType, int order, string sptsType);

//! Get the point locations of the requested type (i.e. Gauss, Lobatto) for the given order
vector<double> getPts1D(string ptsType, int order);

//! Get the Gauss quadrature weights for the Gauss points of the given order [2D]
vector<double> getQptWeights(int order, int nDims);

//! Get the Gauss quadrature weights for the Gauss points of the given order [1D]
vector<double> getQptWeights1D(int order);

//! Get quadrature rule (points & weights) for a tetrahedron for a given order
void getQuadRuleTet(int order, vector<point> &locQpts, vector<double> &weights);

//! Get quadrature rule (points & weights) for a triangle for a given order
void getQuadRuleTri(int order, vector<point> &locQpts, vector<double> &weights);
