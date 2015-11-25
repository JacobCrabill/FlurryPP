/*!
 * \file funcs.hpp
 * \brief Miscellaneous helper functions (header)
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

#include "global.hpp"
#include "input.hpp"

/* --- Linear-Algebra Routines --- */

vector<double> solveCholesky(matrix<double> A, vector<double> b);

matrix<double> solveCholesky(matrix<double> A, matrix<double> &B);

/* ---- Nodal Shape Functions ---- */

//! Shape function for linear or quadratic quad (TODO: Generalize to N-noded quad)
void shape_quad(const point &in_rs, vector<double> &out_shape, int nNodes);
void shape_quad(const point &in_rs, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic quad
void dshape_quad(const point &in_rs, matrix<double> &out_dshape, int nNodes);

//! Second derivative (Hessian) of shape functions for linear or quadratic quad
void ddshape_quad(const point &in_rs, Array<double,3> &out_dshape, int nNodes);

//! Shape function for linear or quadratic hexahedron
void shape_hex(const point &in_rst, vector<double> &out_shape, int nNodes);
void shape_hex(const point &in_rst, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic hexahedron
void dshape_hex(const point &in_rst, matrix<double> &out_dshape, int nNodes);

//! Shape function for linear triangle (TODO: Generalize to N-noded tri)
void shape_tri(const point &in_rs, vector<double> &out_shape);
void shape_tri(const point &in_rs, double* out_shape);

//! Derivative of shape functions for linear triangle
void dshape_tri(point &in_rs, matrix<double> &out_dshape);

//! Shape function for linear tetrahedron
void shape_tet(const point &in_rs, vector<double> &out_shape);
void shape_tet(const point &in_rs, double* out_shape);

//! Derivative of shape functions for linear tetrahedron
void dshape_tet(point &in_rs, matrix<double> &out_dshape);


/* ---- Other ---- */

void getSimplex(int nDims, vector<double> x0, double L, matrix<double> X);

vector<int> getOrder(vector<double> &data);

//! Given points for a cell's face and a point inside the cell, get the outward unit normal
Vec3 getFaceNormalTri(vector<point> &facePts, point &xc);

//! Given points for a cell's face and a point inside the cell, get the outward unit normal
Vec3 getFaceNormalQuad(vector<point> &facePts, point &xc);

//! Given a 2D edge and a point inside the cell, get the outward unit normal
Vec3 getEdgeNormal(vector<point> &edge, point &xc);

//! Given a list of points, get an axis-aligned bounding box defined by its centriod and dimensions
void getBoundingBox(vector<point> &pts, point &minPt, point &maxPt);
void getBoundingBox(matrix<double> &pts, point &minPt, point &maxPt);
void getBoundingBox(double *pts, int nPts, int nDims, point &minPt, point &maxPt);
void getBoundingBox(double *pts, int nPts, int nDims, double *bbox);

vector<double> calcError(const vector<double> &U, const point &pos, input *params);

void calcSolutionFromFlux(matrix<double> &F, vector<double> &U, input *params);

void calcFluxJacobian2D(const vector<double> &U, matrix<double> &dFdU, matrix<double> &dGdU, input *params);
