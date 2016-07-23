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

#include <iomanip>

class geo;

#include "global.hpp"
#include "geo.hpp"
#include "input.hpp"

/* --- Linear-Algebra Routines --- */

vector<double> solveCholesky(matrix<double> A, vector<double> b);

matrix<double> solveCholesky(matrix<double> A, matrix<double> &B);

/* ---- Nodal Shape Functions ---- */

//! Shape function for linear or quadratic quad (TODO: Generalize to N-noded quad)
void shape_quad(const point &in_rs, vector<double> &out_shape, int nNodes);
void shape_quad(const point &in_rs, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic quad
void dshape_quad(const vector<point> loc_pts, Array<double,3> &out_dshape, int nNodes);
void dshape_quad(const point &in_rs, matrix<double> &out_dshape, int nNodes);
void dshape_quad(const point &in_rs, double* out_dshape, int nNodes);

//! Second derivative (Hessian) of shape functions for linear or quadratic quad
void ddshape_quad(const point &in_rs, Array<double,3> &out_dshape, int nNodes);

//! Shape function for linear or quadratic hexahedron
void shape_hex(const point &in_rst, vector<double> &out_shape, int nNodes);
void shape_hex(const point &in_rst, double* out_shape, int nNodes);

//! Derivative of shape functions for linear or quadratic hexahedron
void dshape_hex(const vector<point>& loc_pts, Array<double,3> &out_dshape, int nNodes);
void dshape_hex(const point &in_rst, matrix<double> &out_dshape, int nNodes);
void dshape_hex(const point &in_rst, double* out_dshape, int nNodes);

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

vector<double> calcError(const double * const U, const point &pos, input *params);

void calcFluxJacobian2D(const vector<double> &U, matrix<double> &dFdU, matrix<double> &dGdU, input *params);

//! Given basic grid connectivity for quad mesh, refine by splitting
void refineGridBySplitting2D(matrix<int> &c2v, matrix<int> &c2f, matrix<int> &f2v, vector<point> &xv,
                             vector<int> &parentCell, vector<int> &parentFace);

void refineGrid2D(geo &grid_c, geo &grid_f, int nLevels, int nNodes_c, int shapeOrder_f);

//! Get reference location out_rst of point in_xyz within an element defined by the points in xv
bool getRefLocNewton(double *xv, double *in_xyz, double *out_rst, int nNodes, int nDims);

//! Compute the volume of a high-order quad or hex
double computeVolume(double *xv, int nNodes, int nDims);

/*!
 * Nelder-Mead Minimzation Routine.
 *
 * minFunc should be a normal or lambda function accepting a vector<double>
 * and returning a double.
 */
template<typename Func>
vector<double> NelderMead(const vector<double> &U0, Func minFunc)
{
  // Use the simple Nelder-Meade algorithm to find the reference location which
  // maps to the given physical position

  int nVars = U0.size();
  int nPts = nVars+1;
  vector<std::pair<double,vector<double>>> FX(nPts);

  // Starting location for search
  for (int i=0; i<nPts; i++) {
    FX[i].second = U0;
    if (i>0) {
      FX[i].second[i-1] += .03*FX[i].second[i-1];
    } else {
      for (int j=0; j<nVars; j++) {
        FX[i].second[j] -= .01*FX[i].second[j];
      }
    }
  }

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++)
    FX[i].first = minFunc(FX[i].second);

  std::sort(FX.begin(),FX.end());

  // Use a relative tolerance...?
  double tol = 1e-8;
  int iter = 0;
  while (iter < 200 && FX[0].first>tol) {
    vector<double> Xn = FX[nVars].second;  // Point with the highest value of F
    vector<double> X0(nVars);              // Centroid of all other points
    vector<double> Xr(nVars);              // Reflected point

    // Take centroid of all points besides Xn
    for (int j=0; j<nPts-1; j++)
      for (int k=0; k<nVars; k++)
        X0[k] += FX[j].second[k]/(nPts-1);

    // Reflect Xn around X0
    for (int k=0; k<nVars; k++)
      Xr[k] = X0[k] + (X0[k]-Xn[k]);

    double Fr = minFunc(Xr);

    // Determine what to do with the new point
    if (Fr < FX[nPts-2].first) {
      // We will be keeping this point
      if (Fr < FX[0].first) {
        // This one's good; keep going! Expand from Xr
        vector<double> Xe(nVars);
        for (int i=0; i<nVars; i++)
          Xe[i] = Xr[i] + (X0[i]-Xn[i]);
        double Fe = minFunc(Xe);

        if (Fe < Fr) {
          // This one's even better; use it instead
          FX[nPts-1].first = Fe;
          FX[nPts-1].second = Xe;
        }
        else {
          // Xe/Fe was no better; stick with Fr, Xr
          FX[nPts-1].first = Fr;
          FX[nPts-1].second = Xr;
        }
      }
      else {
        // This one's somewhere in the middle; replace Xn with Xr
        FX[nPts-1].first = Fr;
        FX[nPts-1].second = Xr;
      }
    }
    else {
      // Try reducing the size of the simplex
      vector<double> Xc(nVars);
      for (int i=0; i<nVars; i++)
        Xc[i] = X0[i] - (X0[i]-Xn[i])*.5;
      double Fc = minFunc(Xc);
      if (Fc < FX[nPts-1].first) {
        // Bringing this point in is better; use it
        FX[nPts-1].first = Fc;
        FX[nPts-1].second = Xc;
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        vector<double> X1 = FX[0].second;
        for (int i=1; i<nPts; i++) {
          for (int j=0; j<nVars; j++) {
            FX[i].second[j] = X1[j] + 0.5*(FX[i].second[j]-X1[j]);
          }
          FX[i].first = minFunc(FX[i].second);
        }
      }
    }

    std::sort(FX.begin(),FX.end());

    // Continue to iterate
    iter++;
  }

  return FX[0].second;
}
