/*!
 * \file ele.cpp
 * \brief ele class definition
 *
 * Each elements stores its solution and basic properties, like
 * element type, vertices, and polynomial order
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

#include "ele.hpp"

#include <sstream>

#include "polynomials.hpp"
#include "flux.hpp"
#include "funcs.hpp"

using namespace std;

ele::ele()
{

}

ele::ele(int in_eType, int in_order, int in_ID, vector<point> &in_nodes, geo *in_Geo)
{
  eType = in_eType;
  order = in_order;
  ID = in_ID;
  sID = ID;
  Geo = in_Geo;
}

void ele::initialize(void)
{

}

void ele::setup(input *inParams, solver *inSolver, geo *inGeo, int in_order)
{
  /* --- Basic Stuff --- */
  params = inParams;
  Geo = inGeo;
  Solver = inSolver;

  if (in_order >= 0)
    order = in_order;
  else
    order = Solver->order;

  nDims = params->nDims;
  nFields = params->nFields;

  nSpts = Solver->opers[order].nSpts;
  nFpts = Solver->opers[order].nFpts;

  /* --- Setup all data arrays --- */
  setupArrays();

  /* --- Final Step: calculate physical->reference transforms
   * and store shape basis values for future use --- */
  calcTransforms();

}

void ele::setupArrays(void)
{
  nRKSteps = params->nRKSteps;

  waveSp_fpts.resize(nFpts);

  if (params->scFlag)
    sensor = 0;

  if (params->equation == NAVIER_STOKES && params->calcEntropySensor) {
    S_spts.setup(nSpts,1);
    S_fpts.setup(nFpts,1);
    S_mpts.setup(nMpts,1);
  }

  //tempF.setup(nDims,nFields);
  tempU.assign(nFields,0);
}

/* ---- Solution data access ---- */

double& ele::U_spts(int spt, int field)
{
  return Solver->U_spts(spt, sID, field);
}

double& ele::F_spts(int dim, int spt, int field)
{
  return Solver->F_spts(dim, spt, sID, field);
}

double& ele::dU_spts(int dim, int spt, int field)
{
  return Solver->dU_spts(dim, spt, sID, field);
}

double& ele::dF_spts(int dim_grad, int dim_flux, int spt, int field)
{
  return Solver->dF_spts(dim_grad, dim_flux)(spt, sID, field);
}

double& ele::U_fpts(int fpt, int field)
{
  return Solver->U_fpts(fpt, sID, field);
}

double& ele::F_fpts(int dim, int fpt, int field)
{
  return Solver->F_fpts(dim, fpt, sID, field);
}

double& ele::Fn_fpts(int fpt, int field)
{
  return Solver->Fn_fpts(fpt, sID, field);
}

double& ele::disFn_fpts(int fpt, int field)
{
  return Solver->disFn_fpts(fpt, sID, field);
}

double& ele::dU_fpts(int dim, int fpt, int field)
{
  return Solver->dU_fpts(dim, fpt, sID, field);
}

double& ele::dUc_fpts(int fpt, int field)
{
  return Solver->dUc_fpts(fpt, sID, field);
}

double& ele::divF_spts(int step, int spt, int field)
{
  return Solver->divF_spts[step](spt, sID, field);
}

double& ele::U_mpts(int mpt, int field)
{
  return Solver->U_mpts(mpt, sID, field);
}

double& ele::U_ppts(int ppt, int field)
{
  return Solver->U_ppts(ppt, sID, field);
}

/* ---- Geometry-Variable Access ---- */

double& ele::shape_spts(int spt, int node)
{
  return Solver->shape_spts(spt, node);
}

double& ele::shape_fpts(int fpt, int node)
{
  return Solver->shape_fpts(fpt, node);
}

double& ele::shape_ppts(int ppt, int node)
{
  return Solver->shape_ppts(ppt, node);
}

double& ele::dshape_spts(int spt, int node, int dim)
{
  return Solver->dshape_spts(spt, node, dim);
}

double& ele::dshape_fpts(int fpt, int node, int dim)
{
  return Solver->dshape_fpts(fpt, node, dim);
}

double& ele::detJac_spts(int spt)
{
  return Solver->detJac_spts(spt, sID);
}

double& ele::detJac_fpts(int fpt)
{
  return Solver->detJac_fpts(fpt, sID);
}

double& ele::dA_fpts(int fpt)
{
  return Solver->dA_fpts(fpt, sID);
}

double& ele::Jac_spts(int spt, int dim1, int dim2)
{
  return Solver->Jac_spts(spt, sID, dim1, dim2);
}

double& ele::Jac_fpts(int fpt, int dim1, int dim2)
{
  return Solver->Jac_fpts(fpt, sID, dim1, dim2);
}

double& ele::JGinv_spts(int spt, int dim1, int dim2)
{
  return Solver->JGinv_spts(spt, sID, dim1, dim2);
}

double& ele::JGinv_fpts(int fpt, int dim1, int dim2)
{
  return Solver->JGinv_fpts(fpt, sID, dim1, dim2);
}

double& ele::norm_fpts(int fpt, int dim)
{
  return Solver->norm_fpts(fpt, sID, dim);
}

double& ele::tNorm_fpts(int fpt, int dim)
{
  return Solver->tNorm_fpts(fpt, dim);
}

double& ele::nodes(int npt, int dim)
{
  return Solver->nodes(npt, sID, dim);
}

double& ele::nodesRK(int npt, int dim)
{
  return Solver->nodesRK(npt, sID, dim);
}

double& ele::pos_spts(int spt, int dim)
{
  return Solver->pos_spts(spt, sID, dim);
}

double& ele::pos_fpts(int fpt, int dim)
{
  return Solver->pos_fpts(fpt, sID, dim);
}

double& ele::pos_ppts(int ppt, int dim)
{
  return Solver->pos_ppts(ppt, sID, dim);
}

double& ele::gridVel_spts(int spt, int dim)
{
  return Solver->gridV_spts(spt, sID, dim);
}

double& ele::gridVel_fpts(int fpt, int dim)
{
  return Solver->gridV_fpts(fpt, sID, dim);
}

double& ele::gridVel_nodes(int mpt, int dim)
{
  return Solver->gridV_mpts(mpt, sID, dim);
}

void ele::calcTransforms(bool moving)
{
  /* --- Calculate Transformation at Solution Points --- */
  for (int spt=0; spt<nSpts; spt++)
  {
    for (int dim1=0; dim1<nDims; dim1++)
      for (int dim2=0; dim2<nDims; dim2++)
        Jac_spts(spt,dim1,dim2) = 0.;

    if (!moving)
    {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_spts(spt,dim1,dim2) += dshape_spts(spt,i,dim2)*nodes(i,dim1);
    }
    else
    {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_spts(spt,dim1,dim2) += dshape_spts(spt,i,dim2)*nodesRK(i,dim1);
    }

    if (nDims == 2) {
      // Determinant of transformation matrix
      detJac_spts(spt) = Jac_spts(spt,0,0)*Jac_spts(spt,1,1)-Jac_spts(spt,1,0)*Jac_spts(spt,0,1);
      // Inverse of transformation matrix (times its determinant)
      JGinv_spts(spt,0,0) = Jac_spts(spt,1,1);  JGinv_spts(spt,0,1) =-Jac_spts(spt,0,1);
      JGinv_spts(spt,1,0) =-Jac_spts(spt,1,0);  JGinv_spts(spt,1,1) = Jac_spts(spt,0,0);
    }
    else if (nDims == 3) {
      double xr = Jac_spts(spt,0,0);   double xs = Jac_spts(spt,0,1);   double xt = Jac_spts(spt,0,2);
      double yr = Jac_spts(spt,1,0);   double ys = Jac_spts(spt,1,1);   double yt = Jac_spts(spt,1,2);
      double zr = Jac_spts(spt,2,0);   double zs = Jac_spts(spt,2,1);   double zt = Jac_spts(spt,2,2);
      detJac_spts(spt) = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);

      JGinv_spts(spt,0,0) = ys*zt - yt*zs;  JGinv_spts(spt,0,1) = xt*zs - xs*zt;  JGinv_spts(spt,0,2) = xs*yt - xt*ys;
      JGinv_spts(spt,1,0) = yt*zr - yr*zt;  JGinv_spts(spt,1,1) = xr*zt - xt*zr;  JGinv_spts(spt,1,2) = xt*yr - xr*yt;
      JGinv_spts(spt,2,0) = yr*zs - ys*zr;  JGinv_spts(spt,2,1) = xs*zr - xr*zs;  JGinv_spts(spt,2,2) = xr*ys - xs*yr;
    }
    if (detJac_spts(spt)<0) FatalError("Negative Jacobian at solution points.");
  }

  /* --- Calculate Transformation at Flux Points --- */
  for (int fpt=0; fpt<nFpts; fpt++)
  {
    for (int dim1=0; dim1<nDims; dim1++)
      for (int dim2=0; dim2<nDims; dim2++)
        Jac_fpts(fpt,dim1,dim2) = 0.;

    // Calculate transformation Jacobian matrix - [dx/dr, dx/ds; dy/dr, dy/ds]
    if (!moving)
    {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_fpts(fpt,dim1,dim2) += dshape_fpts(fpt,i,dim2)*nodes(i,dim1);
    }
    else
    {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_fpts(fpt,dim1,dim2) += dshape_fpts(fpt,i,dim2)*nodesRK(i,dim1);
    }

    if (nDims == 2) {
      detJac_fpts(fpt) = Jac_fpts(fpt,0,0)*Jac_fpts(fpt,1,1)-Jac_fpts(fpt,1,0)*Jac_fpts(fpt,0,1);
      // Inverse of transformation matrix (times its determinant)
      JGinv_fpts(fpt,0,0) = Jac_fpts(fpt,1,1);  JGinv_fpts(fpt,0,1) =-Jac_fpts(fpt,0,1);
      JGinv_fpts(fpt,1,0) =-Jac_fpts(fpt,1,0);  JGinv_fpts(fpt,1,1) = Jac_fpts(fpt,0,0);
    }
    else if (nDims == 3) {
      double xr = Jac_fpts(fpt,0,0);   double xs = Jac_fpts(fpt,0,1);   double xt = Jac_fpts(fpt,0,2);
      double yr = Jac_fpts(fpt,1,0);   double ys = Jac_fpts(fpt,1,1);   double yt = Jac_fpts(fpt,1,2);
      double zr = Jac_fpts(fpt,2,0);   double zs = Jac_fpts(fpt,2,1);   double zt = Jac_fpts(fpt,2,2);
      detJac_fpts(fpt) = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);
      // Inverse of transformation matrix (times its determinant)
      JGinv_fpts(fpt,0,0) = ys*zt - yt*zs;  JGinv_fpts(fpt,0,1) = xt*zs - xs*zt;  JGinv_fpts(fpt,0,2) = xs*yt - xt*ys;
      JGinv_fpts(fpt,1,0) = yt*zr - yr*zt;  JGinv_fpts(fpt,1,1) = xr*zt - xt*zr;  JGinv_fpts(fpt,1,2) = xt*yr - xr*yt;
      JGinv_fpts(fpt,2,0) = yr*zs - ys*zr;  JGinv_fpts(fpt,2,1) = xs*zr - xr*zs;  JGinv_fpts(fpt,2,2) = xr*ys - xs*yr;
    }

    /* --- Calculate outward unit normal vector at flux point --- */
    // Transform face normal from reference to physical space [JGinv .dot. tNorm]
    for (int dim1=0; dim1<nDims; dim1++) {
      norm_fpts(fpt,dim1) = 0.;
      for (int dim2=0; dim2<nDims; dim2++) {
        norm_fpts(fpt,dim1) += JGinv_fpts(fpt,dim2,dim1) * tNorm_fpts(fpt,dim2);
      }
    }

    // Store magnitude of face normal (equivalent to face area in finite-volume land)
    dA_fpts(fpt) = 0;
    for (int dim=0; dim<nDims; dim++)
      dA_fpts(fpt) += norm_fpts(fpt,dim)*norm_fpts(fpt,dim);
    dA_fpts(fpt) = sqrt(dA_fpts(fpt));

    // Normalize
    // If we have a collapsed edge, the dA will be 0, so just set the normal to 0
    // (A normal vector at a point doesn't make sense anyways)
    if (std::fabs(dA_fpts(fpt)) < 1e-10) {
      dA_fpts(fpt) = 0.;
      for (int dim=0; dim<nDims; dim++)
        norm_fpts(fpt,dim) = 0;
    }
    else {
      for (int dim=0; dim<nDims; dim++)
        norm_fpts(fpt,dim) /= dA_fpts(fpt);
    }
  }
}

void ele::calcTransforms_point(matrix<double> &jacobian, matrix<double> &JGinv, double &detJac, const point &loc)
{
  // Static transform, or space-time transform
  if (params->motion)
    jacobian.setup(nDims+1,nDims+1);
  else
    jacobian.setup(nDims,nDims);

  jacobian.initializeToZero();

  matrix<double> dshape;
  if (nDims==2)
    dshape_quad(loc, dshape, nNodes);
  else
    dshape_hex(loc, dshape, nNodes);

  if (!params->motion) {
    for (int i=0; i<nNodes; i++)
      for (int dim1=0; dim1<nDims; dim1++)
        for (int dim2=0; dim2<nDims; dim2++)
          jacobian(dim1,dim2) += dshape(i,dim2)*nodes(i,dim1);
  }
  else {
    vector<double> shape;
    shape_quad(loc, shape, nNodes);
    for (int i=0; i<nNodes; i++) {
      for (int dim1=0; dim1<nDims; dim1++) {
        for (int dim2=0; dim2<nDims; dim2++) {
          jacobian(dim1,dim2) += dshape(i,dim2)*nodesRK(i,dim1);
        }
        jacobian(dim1,nDims) += shape[i]*gridVel_nodes(i,dim1);
      }
    }
    jacobian(nDims,nDims) = 1;
  }

  // Determinant of transformation matrix
  detJac = jacobian.det();

  // Inverse of transformation matrix (times its determinant)
  JGinv = jacobian.adjoint();

  if (detJac<0) FatalError("Negative Jacobian at given point.");
}

point ele::calcPos(const point &loc)
{
  getShape(loc,tmpShape);

  point pt;
  if (params->motion == 0) {
    for (int iv=0; iv<nNodes; iv++)
      for (int dim=0; dim<nDims; dim++)
        pt[dim] += tmpShape[iv]*nodes(iv,dim);
  }
  else {
    for (int iv=0; iv<nNodes; iv++)
      for (int dim=0; dim<nDims; dim++)
        pt[dim] += tmpShape[iv]*nodesRK(iv,dim);
  }

  return pt;
}

vector<double> ele::getBoundingBox(void)
{
  vector<double> bbox = {INFINITY,INFINITY,INFINITY,-INFINITY,-INFINITY,-INFINITY};
  if (params->motion == 0) {
    for (uint npt = 0; npt < nNodes; npt++) {
      point pt = point(&nodes(npt,0),nDims);
      for (int dim=0; dim<3; dim++) {
        bbox[dim]   = min(bbox[dim],  pt[dim]);
        bbox[dim+3] = max(bbox[dim+3],pt[dim]);
      }
    }
  }
  else {
    for (uint npt = 0; npt < nNodes; npt++) {
      point pt = point(&nodesRK(npt,0),nDims);
      for (int dim=0; dim<3; dim++) {
        bbox[dim]   = min(bbox[dim],  pt[dim]);
        bbox[dim+3] = max(bbox[dim+3],pt[dim]);
      }
    }
  }

  return bbox;
}

bool ele::getRefLocNewton(point pos, point &loc)
{
  // First, do a quick check to see if the point is even close to being in the element
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  xmin = ymin = zmin =  1e15;
  xmax = ymax = zmax = -1e15;
  double eps = 1e-10;

  auto box = getBoundingBox();
  xmin = box[0];  ymin = box[1];  zmin = box[2];
  xmax = box[3];  ymax = box[4];  zmax = box[5];

  if (pos.x < xmin-eps || pos.y < ymin-eps || pos.z < zmin-eps ||
      pos.x > xmax+eps || pos.y > ymax+eps || pos.z > zmax+eps) {
    // Point does not lie within cell - return an obviously bad ref position
    loc = {99.,99.,99.};
    return false;
  }

  // Use a relative tolerance to handle extreme grids
  double h = min(xmax-xmin,ymax-ymin);
  if (nDims==3) h = min(h,zmax-zmin);

  double tol = 1e-12*h;

  vector<double> shape(nNodes);
  matrix<double> dshape(nNodes,nDims);
  matrix<double> grad(nDims,nDims);

  int iter = 0;
  int iterMax = 20;
  double norm = 1;
  loc = {0, 0, 0};
  while (norm > tol && iter<iterMax) {
    shape_quad(loc,shape,nNodes);
    dshape_quad(loc,dshape,nNodes);

    point dx = pos;
    grad.initializeToZero();
    if (params->motion) {
      for (int n=0; n<nNodes; n++) {
        for (int i=0; i<nDims; i++) {
          for (int j=0; j<nDims; j++) {
            grad(i,j) += nodesRK(n,i)*dshape(n,j);
          }
          dx[i] -= shape[n]*nodesRK(n,i);
        }
      }
    } else {
      for (int n=0; n<nNodes; n++) {
        for (int i=0; i<nDims; i++) {
          for (int j=0; j<nDims; j++) {
            grad(i,j) += nodes(n,i)*dshape(n,j);
          }
          dx[i] -= shape[n]*nodes(n,i);
        }
      }
    }

    double detJ = grad.det();

    auto ginv = grad.adjoint();

    point delta = {0,0,0};
    for (int i=0; i<nDims; i++)
      for (int j=0; j<nDims; j++)
        delta[i] += ginv(i,j)*dx[j]/detJ;

    norm = 0;
    for (int i=0; i<nDims; i++) {
      norm += dx[i]*dx[i];
      loc[i] += delta[i];
      loc[i] = max(min(loc[i],1.),-1.);
    }

    iter++;
    if (iter == iterMax) {
      return false;
    }
  }

  return true;
}

double ele::getDxNelderMead(point refLoc, point physPos)
{
  point pt = calcPos(refLoc);
  Vec3 dx = physPos - pt;

  double norm = dx.norm();

  refLoc.abs();
  for (int i=0; i<nDims; i++) {
    if (refLoc[i]>1) {
      double dxi2 = (refLoc[i]-1.)*(refLoc[i]-1.);
      norm += std::exp(dxi2*dxi2) - 1.;
    }
  }

  return norm;
}

bool ele::getRefLocNelderMead(point pos, point& loc)
{
  // First, do a quick check to see if the point is even close to being in the element
  double xmin, ymin, zmin;
  double xmax, ymax, zmax;
  xmin = ymin = zmin =  1e15;
  xmax = ymax = zmax = -1e15;
  double eps = 1e-10;

  auto box = getBoundingBox();
  xmin = box[0];  ymin = box[1];  zmin = box[2];
  xmax = box[3];  ymax = box[4];  zmax = box[5];

  if (pos.x < xmin-eps || pos.y < ymin-eps || pos.z < zmin-eps ||
      pos.x > xmax+eps || pos.y > ymax+eps || pos.z > zmax+eps) {
    // Point does not lie within cell - return an obviously bad ref position
    loc.x = 99; loc.y = 99; loc.z = 99;
    return false;
  }

  // Use the simple Nelder-Meade algorithm to find the reference location which
  // maps to the given physical position

  int nPts = nDims+1;
  int nVars = nDims;
  vector<std::pair<double,point>> FX(nPts);

  // Starting location for search
  double L = .75;
  if (nDims == 3) {
    FX[0].second = point(-L*.5, -L*.43301, -L*.375);
    FX[1].second = point( L*.5, -L*.43301, -L*.375);
    FX[2].second = point( L*0.,  L*.43301, -L*.375);
    FX[3].second = point( L*0., -L*0.,      L*.375);
  }
  else {
    FX[0].second = point(-L*.5, -L*.43301, 0);
    FX[1].second = point( L*.5, -L*.43301, 0);
    FX[2].second = point( L*0.,  L*.43301, 0);
  }

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++)
    FX[i].first = getDxNelderMead(FX[i].second,pos);

  std::sort(FX.begin(),FX.end());

  // Use a relative tolerance to handle extreme grids
  double h = min(xmax-xmin,ymax-ymin);
  if (nDims==3) h = min(h,zmax-zmin);

  double tol = 1e-10*h;
  int iter = 0;
  while (iter < 300 && FX[0].first>tol) {
    point Xn = FX[nPts-1].second;  // Point with the highest value of F
    point X0;                      // Centroid of all other points
    point Xr;                      // Reflected point

    // Take centroid of all points besides Xn
    for (int j=0; j<nPts-1; j++)
      X0 += FX[j].second/(nPts-1);
    // Reflect Xn around X0
    Xr = X0 + (X0-Xn);

    double Fr = getDxNelderMead(Xr,pos);

    // Determine what to do with the new point
    if (Fr < FX[nPts-2].first) {
      // We will be keeping this point
      if (Fr < FX[0].first) {
        // This one's good; keep going! Expand from Xr
        point Xe = Xr + (X0-Xn);
        double Fe = getDxNelderMead(Xe,pos);

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
      point Xc = X0 - (X0-Xn)*.5;
      double Fc = getDxNelderMead(Xc,pos);
      if (Fc < FX[nPts-1].first) {
        // Bringing this point in is better; use it
        FX[nPts-1].first = Fc;
        FX[nPts-1].second = Xc;
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        point X1 = FX[0].second;
        for (int i=1; i<nPts; i++) {
          for (int j=0; j<nVars; j++) {
            FX[i].second[j] = X1[j] + 0.5*(FX[i].second[j]-X1[j]);
          }
          FX[i].first = getDxNelderMead(FX[i].second,pos);
        }
      }
    }

    std::sort(FX.begin(),FX.end());

    // Continue to iterate
    iter++;
  }

  loc = FX[0].second;

  // Check to see if final location lies within element or not
  eps = 1e-6;
  if (std::abs(loc.x)-eps<=1 && std::abs(loc.y)-eps<=1 && std::abs(loc.z)-eps<=1 && !std::isnan(loc.norm()))
    return true;
  else
    return false;
}

//void ele::calcPosSpts(void)
//{
//  for (int spt=0; spt<nSpts; spt++) {
//    pos_spts(spt,2)ero();
//    for (int iv=0; iv<nNodes; iv++) {
//      for (int dim=0; dim<nDims; dim++) {
//        pos_spts[spt][dim] += shape_spts(spt,iv)*nodes[iv][dim];
//      }
//    }
//  }
//}

//void ele::calcPosFpts(void)
//{
//  for (int fpt=0; fpt<nFpts; fpt++) {
//    pos_fpts[fpt].zero();
//    for (int iv=0; iv<nNodes; iv++) {
//      for (int dim=0; dim<nDims; dim++) {
//        pos_fpts[fpt][dim] += shape_fpts(fpt,iv)*nodes[iv][dim];
//      }
//    }
//  }
//}

//void ele::updatePosSpts(void)
//{
//  for (int spt=0; spt<nSpts; spt++) {
//    pos_spts(spt,2)ero();
//    for (int iv=0; iv<nNodes; iv++) {
//      for (int dim=0; dim<nDims; dim++) {
//        pos_spts[spt][dim] += shape_spts(spt,iv)*nodesRK[iv][dim];
//      }
//    }
//  }
//}

//void ele::updatePosFpts(void)
//{
//  for (int fpt=0; fpt<nFpts; fpt++) {
//    pos_fpts[fpt].zero();
//    for (int iv=0; iv<nNodes; iv++) {
//      for (int dim=0; dim<nDims; dim++) {
//        pos_fpts[fpt][dim] += shape_fpts(fpt,iv)*nodesRK[iv][dim];
//      }
//    }
//  }
//}


void ele::setInitialCondition()
{
  if (params->equation == NAVIER_STOKES) {
    double rho, vx, vy, vz, p;
    double gamma = params->gamma;

    if (params->icType == 0) {
      /* --- Uniform "Freestream" solution --- */
      rho = params->rhoIC;
      vx = params->vxIC;
      vy = params->vyIC;
      if (nDims == 3)
        vz = params->vzIC;
      else
        vz = 0;

      p = params->pIC;
      for (int spt=0; spt<nSpts; spt++) {
        U_spts(spt,0) = rho;
        U_spts(spt,1) = rho * vx;
        U_spts(spt,2) = rho * vy;
        if (nDims == 3) U_spts(spt,3) = rho * vz;
        U_spts(spt,nDims+1) = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy + vz*vz));
      }
    }
    else if (params->icType == 1) {
      /* --- Isentropic Vortex of strength eps centered at (0,0) --- */
      double eps = 5.0;
      for (int spt=0; spt<nSpts; spt++) {
        double x = pos_spts(spt,0);
        double y = pos_spts(spt,1);

        double f = 1.0 - (x*x + y*y);

        // Limiting rho to 1e-3 to avoid negative density/pressure issues
        rho = max(pow(1. - eps*eps*(gamma-1.)/(8.*gamma*pi*pi)*exp(f), 1.0/(gamma-1.0) + 1e-5), 1e-3);
        vx = 1. - eps*y / (2.*pi) * exp(f/2.);
        vy = 1. + eps*x / (2.*pi) * exp(f/2.);
        p = pow(rho,gamma);

        U_spts(spt,0) = rho;
        U_spts(spt,1) = rho * vx;
        U_spts(spt,2) = rho * vy;
        if (nDims == 3) U_spts(spt,3) = 0.;
        U_spts(spt,nDims+1) = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
    else if (params->icType == 2) {
      /* --- Isentropic Vortex of strength eps centered at (0,0) (Liang version) --- */
      double eps = 1.0;  // See paper by Liang and Miyaji, CPR Deforming Domains
      double rc  = 1.0;
      double Minf = .3;
      double Uinf = 1;
      double rhoInf = 1;
      double theta = atan(0.5);
      double Pinf = pow(Minf,-2)/gamma;

      double eM = (eps*Minf)*(eps*Minf);

      for (int spt=0; spt<nSpts; spt++) {
        double x = pos_spts(spt,0);
        double y = pos_spts(spt,1);

        double f = -(x*x + y*y) / (rc*rc);

        vx = Uinf*(cos(theta) - y*eps/rc * exp(f/2.));
        vy = Uinf*(sin(theta) + x*eps/rc * exp(f/2.));
        rho = rhoInf*pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));
        p   = Pinf  *pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));

        U_spts(spt,0) = rho;
        U_spts(spt,1) = rho * vx;
        U_spts(spt,2) = rho * vy;
        if (nDims == 3) U_spts(spt,3) = 0.;
        U_spts(spt,nDims+1) = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    if (params->icType == 0) {
      /* --- Simple Gaussian bump centered at (0,0) --- */
      for (int spt=0; spt<nSpts; spt++) {
        point pt = point(&pos_spts(spt,0),nDims);
        double r2 = pt*pt;
        U_spts(spt,0) = exp(-r2);
      }
    }
    else if (params->icType == 1) {
      /* --- Test case: sin(x) --- */
      for (int spt=0; spt<nSpts; spt++) {
        U_spts(spt,0) = 1. + sin(2.*pi*(pos_spts(spt,0)+5)/10.);
      }
    }
    else if (params->icType == 2) {
      /* --- Test case for debugging - cos(x)*cos(y)*cos(z) over domain --- */
      for (int spt=0; spt<nSpts; spt++)
        U_spts(spt,0) = cos(2*pi*pos_spts(spt,0)/6.)*cos(2*pi*pos_spts(spt,1)/6.)*cos(2*pi*pos_spts(spt,2)/6.);
    }
  }
}

matrix<double> ele::calcError(void)
{

  matrix<double> err(nSpts,nFields);

  if (!params->testCase) {
    for (uint spt = 0; spt < nSpts; spt++)
      for (uint k = 0; k < nFields; k++)
        err(spt,k) = U_spts(spt,k);
    return err;
  }

  if (params->equation == NAVIER_STOKES) {
    double gamma = params->gamma;

    if (params->icType == 0) {
      /* --- Uniform "Freestream" solution --- */
      double rho = params->rhoIC;
      double vx = params->vxIC;
      double vy = params->vyIC;
      double vz = 0.;
      if (nDims == 3)
        vz = params->vzIC;

      double p = params->pIC;
      for (int spt=0; spt<nSpts; spt++) {
        err(spt,0) = rho;
        err(spt,1) = rho * vx;
        err(spt,2) = rho * vy;
        if (nDims == 3) err(spt,3) = rho * vz;
        err(spt,nDims+1) = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy + vz*vz));
      }
    }
    else if (params->icType == 1) {
      /* --- Isentropic Vortex of strength eps centered at (0,0) --- */
      double eps = 5.0;

      double xmin, xmax, ymin, ymax;
      if (params->meshType == CREATE_MESH) {
        xmin = params->xmin;
        xmax = params->xmax;
        ymin = params->ymin;
        ymax = params->ymax;
      } else {
        // Assuming a 'standard' mesh for the test case
        xmin = -5;  xmax = 5;
        ymin = -5;  ymax = 5;
      }

      for (int spt=0; spt<nSpts; spt++) {
        double x = fmod( (pos_spts(spt,0) - params->time), (xmax-xmin) );
        double y = fmod( (pos_spts(spt,1) - params->time), (ymax-ymin) );
        if (x > xmax) x -= (xmax-xmin);
        if (y > ymax) y -= (ymax-ymin);
        if (x < xmin) x += (xmax-xmin);
        if (y < ymin) y += (ymax-ymin);

        double f = 1.0 - (x*x + y*y);

        // Limiting rho to 1e-3 to avoid negative density/pressure issues
        double rho = max(pow(1. - eps*eps*(gamma-1.)/(8.*gamma*pi*pi)*exp(f), 1.0/(gamma-1.0) + 1e-5), 1e-3);
        double vx = 1. - eps*y / (2.*pi) * exp(f/2.);
        double vy = 1. + eps*x / (2.*pi) * exp(f/2.);
        double p = pow(rho,gamma);

        err(spt,0) = rho;
        err(spt,1) = rho * vx;
        err(spt,2) = rho * vy;
        if (nDims == 3) err(spt,3) = 0.;
        err(spt,nDims+1) = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
    else if (params->icType == 2) {
      /* --- Isentropic Vortex of strength eps centered at (0,0) (Liang version) --- */
      double eps = 1.0;  // See paper by Liang and Miyaji, CPR Deforming Domains
      double rc  = 1.0;
      double Minf = .3;
      double Uinf = 1;
      double rhoInf = 1;
      double theta = atan(0.5);
      double Pinf = pow(Minf,-2)/gamma;

      double eM = (eps*Minf)*(eps*Minf);

      double xmin, xmax, ymin, ymax;
      if (params->meshType == CREATE_MESH) {
        xmin = params->xmin;
        xmax = params->xmax;
        ymin = params->ymin;
        ymax = params->ymax;
      } else {
        // Assuming a 'standard' mesh for the test case
        xmin = -5;  xmax = 5;
        ymin = -5;  ymax = 5;
      }

      for (int spt=0; spt<nSpts; spt++) {
        double x = fmod( (pos_spts(spt,0) - Uinf*cos(theta)*params->time), (xmax-xmin) );
        double y = fmod( (pos_spts(spt,1) - Uinf*sin(theta)*params->time), (ymax-ymin) );
        if (x > xmax) x -= (xmax-xmin);
        if (y > ymax) y -= (ymax-ymin);
        if (x < xmin) x += (xmax-xmin);
        if (y < ymin) y += (ymax-ymin);

        double f = -(x*x + y*y) / (rc*rc);

        double vx = Uinf*(cos(theta) - y*eps/rc * exp(f/2.));
        double vy = Uinf*(sin(theta) + x*eps/rc * exp(f/2.));
        double rho = rhoInf*pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));
        double p   = Pinf  *pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));

        err(spt,0) = rho;
        err(spt,1) = rho * vx;
        err(spt,2) = rho * vy;
        if (nDims == 3) err(spt,3) = 0.;
        err(spt,nDims+1) =  p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
      }
    }
  }
  else if (params->equation == ADVECTION_DIFFUSION) {
    double xmin, xmax, ymin, ymax;
    if (params->meshType == CREATE_MESH) {
      xmin = params->xmin;
      xmax = params->xmax;
      ymin = params->ymin;
      ymax = params->ymax;
    } else {
      // Assuming a 'standard' mesh for the test case
      xmin = -5;  xmax = 5;
      ymin = -5;  ymax = 5;
    }

    if (params->icType == 0) {
      /* --- Simple Gaussian bump centered at (0,0) --- */
      for (int spt=0; spt<nSpts; spt++) {
        double x = std::abs(fmod( (pos_spts(spt,0) - params->time), (xmax-xmin) ));
        double y = std::abs(fmod( (pos_spts(spt,1) - params->time), (ymax-ymin) ));
        if (x > xmax) x -= (xmax-xmin);
        if (y > ymax) y -= (ymax-ymin);

        double r2 = x*x + y*y;
        err(spt,0) = exp(-r2);
      }
    }
    else if (params->icType == 1) {
      /* --- Test case equivalent to 1D test - Advection of sine wave --- */
      for (int spt=0; spt<nSpts; spt++)
        err(spt,0) = 1 + sin(2*pi*(pos_spts(spt,0)+5-params->time)/10.);
    }
    else if (params->icType == 2) {
      /* --- Test case for debugging - cos(x)*cos(y)*cos(z) over domain --- */
      for (int spt=0; spt<nSpts; spt++)
        err(spt,0) = cos(2*pi*pos_spts(spt,0)/6.)*cos(2*pi*pos_spts(spt,1)/6.)*cos(2*pi*pos_spts(spt,2)/6.);
    }
  }

  for (int spt=0; spt<nSpts; spt++)
    for (int j=0; j<nFields; j++)
      err(spt,j) = U_spts(spt,j) - err(spt,j);

  if (params->errorNorm == 1)
    for (auto &val:err.data) val = abs(val); // L1 norm
  else if (params->errorNorm == 2)
    for (auto &val:err.data) val *= val;     // L2 norm

  return err;
}

void ele::getShape(point loc, vector<double> &shape)
{
  if (eType == TRI) {
    shape_tri(loc, shape);
  }
  else if (eType == QUAD) {
    shape_quad(loc, shape, nNodes);
  }
  else if (eType == HEX) {
    shape_hex(loc, shape, nNodes);
  }
  else {
    FatalError("Element Type Not Supported.");
  }
}

void ele::calcViscousFlux_spts()
{
  for (int spt=0; spt<nSpts; spt++) {

    // TEMP HACK (inefficient, but will work fine)
    matrix<double> tempDU(nDims,nFields);
    for (int dim=0; dim<nDims; dim++) {
      for (int k=0; k<nFields; k++) {
        tempDU(dim,k) = dU_spts(dim,spt,k);
      }
    }

    if (params->equation == NAVIER_STOKES)
      viscousFlux(&Solver->U_spts(ID,spt), tempDU, tempF, params);
    else if (params->equation == ADVECTION_DIFFUSION)
      viscousFluxAD(tempDU, tempF, params);

    if (params->motion) {
      /* --- Don't transform yet; that will be handled later --- */
      for (int dim=0; dim<nDims; dim++) {
        for (int k=0; k<nFields; k++) {
          F_spts(dim,spt,k) += tempF[dim][k];
        }
      }
    }
    else {
      /* --- Transform back to reference domain --- */
      for (int k=0; k<nFields; k++) {
        for (int dim=0; dim<nDims; dim++) {
          for (int j=0; j<nDims; j++) {
            F_spts(dim,spt,k) += JGinv_spts(spt,dim,j)*tempF[j][k];
          }
        }
      }
    }

  }
}

vector<matrix<double>> ele::transformFlux_physToRef(void)
{
  vector<matrix<double>> outF(nDims);
  for (auto &FD:outF) {
    FD.setup(nSpts,nFields);
    FD.initializeToZero();
  }

  if (params->motion) {
    // Use space-time transformation
    for (int spt=0; spt<nSpts; spt++) {
      matrix<double> jacobian(nDims+1,nDims+1);
      jacobian(nDims,nDims) = 1;
      for (int dim1=0; dim1<nDims; dim1++) {
        jacobian(dim1,nDims) = gridVel_spts(spt,dim1);
        for (int dim2=0; dim2<nDims; dim2++) {
          jacobian(dim1,dim2) = Jac_spts(spt,dim1,dim2);
        }
      }

      auto S = jacobian.adjoint();

      for (int k=0; k<nFields; k++) {
        for (int dim1=0; dim1<nDims; dim1++) {
          outF[dim1](spt,k) = U_spts(spt,k) * S(dim1,nDims);
          for (int dim2=0; dim2<nDims; dim2++) {
            outF[dim1](spt,k) += S(dim1,dim2) * F_spts(dim2,spt,k);
          }
        }
      }
    }
  } else {
    // Standard static transformation
    for (int spt=0; spt<nSpts; spt++) {
      for (int dim1=0; dim1<nDims; dim1++) {
        for (int k=0; k<nFields; k++) {
          outF[dim1](spt,k) = 0.;
          for (int dim2=0; dim2<nDims; dim2++) {
            outF[dim1](spt,k) += JGinv_spts(spt,dim1,dim2)*F_spts(dim2,spt,k);
          }
        }
      }
    }
  }

  return outF;
}

vector<matrix<double>> ele::transformFlux_refToPhys(void)
{
  vector<matrix<double>> outF(nDims);

  for (int spt=0; spt<nSpts; spt++) {
    for (int dim1=0; dim1<nDims; dim1++) {
      for (int k=0; k<nFields; k++) {
        outF[dim1](spt,k) = 0.;
        for (int dim2=0; dim2<nDims; dim2++) {
          outF[dim1](spt,k) += Jac_spts(spt,dim1,dim2) * F_spts(dim2,spt,k) / detJac_spts(spt);
        }
      }
    }
  }

  return outF;
}

vector<matrix<double>> ele::transformGradU_physToRef(void)
{
  vector<matrix<double>> outDU(nDims);
  for (auto &DU:outDU) {
    DU.setup(nSpts,nFields);
    DU.initializeToZero();
  }

  if (nDims == 2) {
    for (int spt=0; spt<nSpts; spt++) {
      for (int k=0; k<nFields; k++) {
        outDU[0](spt,k) =  dU_spts(0,spt,k)*Jac_spts(spt,1,1) - dU_spts(1,spt,k)*Jac_spts(spt,0,1);
        outDU[1](spt,k) = -dU_spts(0,spt,k)*Jac_spts(spt,1,0) + dU_spts(1,spt,k)*Jac_spts(spt,0,0);
      }
    }
  }

  return outDU;
}

void ele::transformGradF_spts(int step)
{
  // NOTE: The 1st dim of dF is the derivative, and the 2nd is the flux direction

  if (nDims == 2) {
    for (int spt=0; spt<nSpts; spt++) {
      double A = gridVel_spts(spt,1)*Jac_spts(spt,0,1) - gridVel_spts(spt,0)*Jac_spts(spt,1,1);
      double B = gridVel_spts(spt,0)*Jac_spts(spt,1,0) - gridVel_spts(spt,1)*Jac_spts(spt,0,0);
      for (int k=0; k<nFields; k++) {
        dF_spts(0,0,spt,k) =  dF_spts(0,0,spt,k)*Jac_spts(spt,1,1) - dF_spts(0,1,spt,k)*Jac_spts(spt,0,1) + dU_spts(0,spt,k)*A;
        dF_spts(1,1,spt,k) = -dF_spts(1,0,spt,k)*Jac_spts(spt,1,0) + dF_spts(1,1,spt,k)*Jac_spts(spt,0,0) + dU_spts(1,spt,k)*B;
        divF_spts(step,spt,k) = dF_spts(0,0,spt,k)+dF_spts(1,1,spt,k);
      }
    }
  } else {
    for (int spt=0; spt<nSpts; spt++) {
      // Build the full 4D (space+time) Jacobian matrix & its adjoint
      matrix<double> Jacobian(4,4);
      Jacobian(3,3) = 1;
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++)
          Jacobian(i,j) = Jac_spts(spt,i,j);
        Jacobian(i,3) = gridVel_spts(spt,i);
      }
      matrix<double> S = Jacobian.adjoint();

      for (int dim1=0; dim1<3; dim1++)
        for (int dim2=0; dim2<3; dim2++)
          for (int k=0; k<nFields; k++)
            divF_spts(step,spt,k) += dF_spts(dim2,dim1,spt,k)*S(dim2,dim1);

      for (int dim=0; dim<3; dim++)
        for (int k=0; k<nFields; k++)
          divF_spts(step,spt,k) += dU_spts(dim,spt,k)*S(dim,3);
    }
  }
}

void ele::calcDeltaUc(void)
{
//  for (int fpt=0; fpt<nFpts; fpt++) {
//    for (int k=0; k<nFields; k++) {
//      dUc_fpts(fpt,k) = Uc_fpts(fpt,k) - U_fpts(fpt,k);
//    }
//  }
}

void ele::calcEntropyErr_spts(void)
{
  for (int spt=0; spt<nSpts; spt++) {
    auto v = getEntropyVars(spt);
    S_spts(spt) = 0;
    for (int k=0; k<nFields; k++) {
      S_spts(spt) += v[k]*divF_spts(0,spt,k);
    }
    S_spts(spt) /= detJac_spts(spt);
  }
}

vector<double> ele::getEntropyVars(int spt)
{
  vector<double> v(nFields);
  double gamma = params->gamma;

  auto phi = getPrimitives(spt);

  if (nDims == 2) {
    double S = log(phi[3]) - gamma*log(phi[0]); // ln(p) - gamma ln(rho)
    double Vmag2 = phi[1]*phi[1] + phi[2]*phi[2];

    v[0] = (gamma-S)/(gamma-1) - 0.5*phi[0]*Vmag2/phi[3];
    v[1] = phi[0]*phi[1]/phi[3];
    v[2] = phi[0]*phi[2]/phi[3];
    v[3] = -phi[0]/phi[3];
  }
  else {
    double S = log(phi[4]) - gamma*log(phi[0]); // ln(p) - gamma ln(rho)
    double Vmag2 = phi[1]*phi[1] + phi[2]*phi[2] + phi[3]*phi[3];

    v[0] = (gamma-S)/(gamma-1) - 0.5*phi[0]*Vmag2/phi[4];
    v[1] = phi[0]*phi[1]/phi[4];
    v[2] = phi[0]*phi[2]/phi[4];
    v[3] = phi[0]*phi[3]/phi[4];
    v[4] = -phi[0]/phi[4];
  }

  return v;
}

void ele::calcWaveSpFpts(void)
{
  if (params->equation == ADVECTION_DIFFUSION) {
    for (int fpt=0; fpt<nFpts; fpt++) {
      double u = params->advectVx;
      double v = params->advectVy;
      double w = 0.;
      if (nDims == 3) w = params->advectVz;

      if (params->motion) {
        u -= gridVel_fpts(fpt,0);
        v -= gridVel_fpts(fpt,1);
        if (nDims == 3)
          w -= gridVel_fpts(fpt,2);
      }

      double csq = u*u + v*v + w*w;

      waveSp_fpts[fpt] = sqrt(csq) / dA_fpts(fpt);
    }
  }
  else if (params->equation == NAVIER_STOKES) {
    for (int fpt=0; fpt<nFpts; fpt++) {
      double rho = U_fpts(fpt,0);
      double u = U_fpts(fpt,1)/rho;
      double v = U_fpts(fpt,2)/rho;
      double w = 0;
      if (nDims == 3) w = U_fpts(fpt,3)/rho;

      double rhoVSq = rho*(u*u+v*v+w*w);
      double p = (params->gamma-1)*(U_fpts(fpt,nDims+1) - 0.5*rhoVSq);

      double vN = u*norm_fpts(fpt,0) + v*norm_fpts(fpt,1);
      if (nDims == 3) vN += w*norm_fpts(fpt,2);

      double vgN = 0;
      if (params->motion) {
        vgN = gridVel_fpts(fpt,0)*norm_fpts(fpt,0) + gridVel_fpts(fpt,1)*norm_fpts(fpt,1);
        if (nDims == 3) vgN += gridVel_fpts(fpt,2)*norm_fpts(fpt,2);
      }

      double csq = std::max(params->gamma*p/rho,0.0);
      waveSp_fpts[fpt] = (std::abs(vN - vgN) + std::sqrt(csq)) / dA_fpts(fpt);
    }
  }
}

double ele::calcDt(void)
{
  double waveSp = 0;
  for (int fpt=0; fpt<nFpts; fpt++)
    if (dA_fpts(fpt) > 0) // ignore collapsed edges
      waveSp = max(waveSp,waveSp_fpts[fpt]);

  dt = (params->CFL) * getCFLLimit(order) * (2.0 / (waveSp+1.e-10));
  return dt;
}

vector<double> ele::getPrimitives(uint spt)
{
  vector<double> V(nFields);

  if (params->equation == ADVECTION_DIFFUSION) {
    V[0] = U_spts(spt,0);
  }
  else if (params->equation == NAVIER_STOKES) {
    V[0] = U_spts(spt,0);
    V[1] = U_spts(spt,1)/V[0];
    V[2] = U_spts(spt,2)/V[0];
    double vMagSq = V[1]*V[1]+V[2]*V[2];
    if (nDims == 3) {
      V[3] = U_spts(spt,3)/V[0];
      vMagSq += V[3]*V[3];
    }
    V[nDims+1] = (params->gamma-1)*(U_spts(spt,nDims+1) - 0.5*V[0]*vMagSq);
  }

  return V;
}

vector<double> ele::getPrimitivesFpt(uint fpt)
{
  vector<double> V(nFields);

  if (params->equation == ADVECTION_DIFFUSION) {
    V[0] = U_fpts(fpt,0);
  }
  else if (params->equation == NAVIER_STOKES) {
    V[0] = U_fpts(fpt,0);
    V[1] = U_fpts(fpt,1)/V[0];
    V[2] = U_fpts(fpt,2)/V[0];
    double vMagSq = V[1]*V[1]+V[2]*V[2];
    if (nDims == 3) {
      V[3] = U_fpts(fpt,3)/V[0];
      vMagSq += V[3]*V[3];
    }
    V[nDims+1] = (params->gamma-1)*(U_fpts(fpt,nDims+1) - 0.5*V[0]*vMagSq);
  }

  return V;
}

vector<double> ele::getPrimitivesMpt(uint mpt)
{
  vector<double> V(nFields);

  if (params->equation == ADVECTION_DIFFUSION) {
    V[0] = U_mpts(mpt,0);
  }
  else if (params->equation == NAVIER_STOKES) {
    V[0] = U_mpts(mpt,0);
    V[1] = U_mpts(mpt,1)/V[0];
    V[2] = U_mpts(mpt,2)/V[0];
    double vMagSq = V[1]*V[1]+V[2]*V[2];
    if (nDims == 3) {
      V[3] = U_mpts(mpt,3)/V[0];
      vMagSq += V[3]*V[3];
    }
    V[nDims+1] = (params->gamma-1)*(U_mpts(mpt,nDims+1) - 0.5*V[0]*vMagSq);
  }

  return V;
}

void ele::getPrimitivesPlot(matrix<double> &V)
{
  V.setup(Solver->nPpts, nFields);

  for (uint ppt = 0; ppt < Solver->nPpts; ppt++)
    for (uint k = 0; k < nFields; k++)
      V(ppt,k) = U_ppts(ppt, k);

  if (params->equation == NAVIER_STOKES) {
    // Overwriting V, so be careful of order!
    for (uint i = 0; i < V.getDim0(); i++) {
      double u = V(i,1) / V(i,0);
      double v = V(i,2) / V(i,0);
      double w = 0.;
      if (nDims == 3) w = V(i,3) / V(i,0);
      double vSq = u*u + v*v + w*w;
      V(i,nDims+1) = (params->gamma-1)*(V(i,nDims+1) - 0.5*V(i,0)*vSq);
      V(i,1) = u;
      V(i,2) = v;
      if (nDims == 3) V(i,3) = w;
    }
  }
}

void ele::getGridVelPlot(matrix<double> &GV)
{
  if (eType == QUAD) {
    GV.setup(nSpts+nFpts+nMpts,nDims);

    // Get solution at corner points
    for (int dim=0; dim<nDims; dim++) {
      GV(0,dim)                     = gridVel_nodes(0,dim);
      GV(order+2,dim)               = gridVel_nodes(1,dim);
      GV((order+3)*(order+3)-1,dim) = gridVel_nodes(2,dim);
      GV((order+3)*(order+2),dim)   = gridVel_nodes(3,dim);
    }

    // Get solution at flux points
    for (int i=0; i<order+1; i++) {
      for (int dim=0; dim<nDims; dim++) {
        GV(i+1,dim)                     = gridVel_fpts(i,dim);               // Bottom
        GV((i+1)*(order+3),dim)         = gridVel_fpts(nFpts-i-1,dim);       // Left
        GV((i+2)*(order+3)-1,dim)       = gridVel_fpts(order+1+i,dim);       // Right
        GV((order+3)*(order+2)+i+1,dim) = gridVel_fpts(3*(order+1)-i-1,dim); // Top
      }
    }

    // Get solution at solution points
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        int id = (i+1)*(order+3)+j+1;
        for (int dim=0; dim<nDims; dim++) {
          GV(id,dim) = gridVel_spts(j+i*(order+1),dim);
        }
      }
    }
  }
  else if (eType == HEX) {
    int nPts1D = order+3;
    int P22 = nPts1D*nPts1D;
    int nv = 8, ne = 12;

    GV.setup(nPts1D*nPts1D*nPts1D,nFields);

    for (int dim=0; dim<nDims; dim++) {
      // Get solution at corner points
      GV(0,dim)                = gridVel_nodes(0,dim);
      GV(order+2,dim)          = gridVel_nodes(1,dim);
      GV(P22-1,dim)            = gridVel_nodes(2,dim);
      GV(nPts1D*(order+2),dim) = gridVel_nodes(3,dim);

      int base = (order+2)*P22;
      GV(base,dim)                    = gridVel_nodes(4,dim);
      GV(base + order+2,dim)          = gridVel_nodes(5,dim);
      GV(base + P22-1,dim)            = gridVel_nodes(6,dim);
      GV(base + nPts1D*(order+2),dim) = gridVel_nodes(7,dim);

      // Get solution at edge points
      for (int i=0; i<order+1; i++) {
        /* --- Bottom Edges --- */
        // edge 0-1
        GV(i+1,dim) = gridVel_nodes(nv+i*ne+0,dim);

        // edge 0-3
        GV(nPts1D*(i+1),dim) = gridVel_nodes(nv+(order-i)*ne+3,dim);

        // edge 1-2
        GV(nPts1D*(i+2)-1,dim) = gridVel_nodes(nv+i*ne+1,dim);

        // edge 3-2
        GV(nPts1D*(order+2)+i+1,dim) = gridVel_nodes(nv+(order-i)*ne+2,dim);

        /* --- Top Edges --- */
        base = P22*(order+2);
        // edge 4-5
        GV(base + i+1,dim) = gridVel_nodes(nv+i*ne+4,dim);

        // edge 4-7
        GV(base + nPts1D*(i+1),dim) = gridVel_nodes(nv+(order-i)*ne+7,dim);

        // edge 5-6
        GV(base + nPts1D*(i+2)-1,dim) = gridVel_nodes(nv+i*ne+5,dim);

        // edge 7-6
        GV(base + nPts1D*(order+2)+i+1,dim) = gridVel_nodes(nv+(order-i)*ne+6,dim);

        /* --- Mid [Vertical] Egdes --- */
        base = (i+1)*P22;
        // edge 0-4
        GV(base,dim) = gridVel_nodes(nv+i*ne+8,dim);

        // edge 1-5
        GV(base+(order+2),dim) = gridVel_nodes(nv+i*ne+9,dim);

        int base2 = nPts1D*(order+2);
        // edge 3-7
        GV(base+base2,dim) = gridVel_nodes(nv+i*ne+11,dim);

        // edge 2-6
        GV(base+base2+order+2,dim) = gridVel_nodes(nv+i*ne+10,dim);
      }

      // Get solution at flux points
      int P12 = (order+1)*(order+1);
      for (int i=0; i<order+1; i++) {
        for (int j=0; j<order+1; j++) {
          int ind1 = i + j*(order+1);
          int ind2 = order - i + (order+1)*j;
          // Bottom Face
          GV(nPts1D*(j+1)+i+1,dim) = gridVel_fpts(ind1,dim);

          // Top Face
          GV(P22*(order+2)+(j+1)*nPts1D+i+1,dim) = gridVel_fpts(P12+ind2,dim);

          // Left Face
          GV(P22*(j+1)+nPts1D*(i+1),dim) = gridVel_fpts(2*P12+ind1,dim);

          // Right Face
          GV(P22*(j+1)+nPts1D*(i+1)+order+2,dim) = gridVel_fpts(3*P12+ind2,dim);

          // Front Face
          GV(P22*(j+1)+i+1,dim) = gridVel_fpts(4*P12+ind2,dim);

          // Back Face
          GV(P22*(j+2)+i+1-nPts1D,dim) = gridVel_fpts(5*P12+ind1,dim);
        }
      }

      // Get solution at solution points
      for (int k=0; k<order+1; k++) {
        for (int j=0; j<order+1; j++) {
          for (int i=0; i<order+1; i++) {
            GV(i+1+nPts1D*(j+1)+(k+1)*P22,dim) = gridVel_spts(i+(order+1)*(j+(order+1)*k),dim);
          }
        }
      }
    }
  }
}

void ele::getEntropyErrPlot(matrix<double> &S)
{
  if (nDims == 3) FatalError("Entropy-error calculation not yet supported for 3D cases.");

  S.setup(nSpts+nFpts+nMpts,1);

  // Get solution at corner points
  S(0)                     = S_mpts(0);
  S(order+2)               = S_mpts(1);
  S((order+3)*(order+3)-1) = S_mpts(2);
  S((order+3)*(order+2))   = S_mpts(3);

  // Get solution at flux points
  for (int i=0; i<order+1; i++) {
      S(i+1)                     = S_fpts(i);               // Bottom
      S((i+1)*(order+3))         = S_fpts(nFpts-i-1);       // Left
      S((i+2)*(order+3)-1)       = S_fpts(order+1+i);       // Right
      S((order+3)*(order+2)+i+1) = S_fpts(3*(order+1)-i-1); // Top
  }

  // Get solution at solution points
  for (int i=0; i<order+1; i++) {
    for (int j=0; j<order+1; j++) {
      int id = (i+1)*(order+3)+j+1;
      S(id) = S_spts(j+i*(order+1));
    }
  }
}

bool ele::checkDensity()
{
  /* --- Fisrt, check if density is negative and squeeze if needed --- */
  bool negRho = false;
  double minRho = 1e15;
  double tol = 1e-10;   // Tolerance for squeezing

  for (int spt=0; spt<nSpts; spt++) {
    if (U_spts(spt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_spts(spt,0));
      break;
    }
  }

  for (int fpt=0; fpt<nFpts; fpt++) {
    if (U_fpts(fpt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_fpts(fpt,0));
      break;
    }
  }

  // --- Do the squeezing on density (if needed) ---
  if (negRho) {
    double eps = abs(Uavg[0] - tol)/(Uavg[0] - minRho);
    for (int spt=0; spt<nSpts; spt++) {
      U_spts(spt,0) = (1-eps)*Uavg[0] + eps*U_spts(spt,0);
    }

    for (int fpt=0; fpt<nFpts; fpt++) {
      U_fpts(fpt,0) = (1-eps)*Uavg[0] + eps*U_fpts(fpt,0);
    }
  }

  return negRho;
}

void ele::checkEntropy()
{
  /* --- Fisrt, check if density is negative and squeeze if needed --- */
  bool negRho = false;
  double minRho = 1e15;
  double tol = 1e-10;   // Tolerance for squeezing

  for (int spt=0; spt<nSpts; spt++) {
    if (U_spts(spt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_spts(spt,0));
    }
  }

  for (int fpt=0; fpt<nFpts; fpt++) {
    if (U_fpts(fpt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_fpts(fpt,0));
    }
  }

  // --- Do the squeezing on density (if needed) ---
  if (negRho) {
    double eps = abs(Uavg[0] - tol)/(Uavg[0] - minRho);
    for (int spt=0; spt<nSpts; spt++) {
      U_spts(spt,0) = (1-eps)*Uavg[0] + eps*U_spts(spt,0);
    }

    for (int fpt=0; fpt<nFpts; fpt++) {
      U_fpts(fpt,0) = (1-eps)*Uavg[0] + eps*U_fpts(fpt,0);
    }
  }

  /* --- Next, check for entropy loss and correct if needed --- */

  double minTau = 1e15; // Entropy-bounding value
  for (int spt=0; spt<nSpts; spt++) {
    auto phi = getPrimitives(spt);
    double rho = phi[0];
    double p = phi[nDims+1];

    // Get minimum 'tau' value
    minTau = std::min(minTau, p - params->exps0*std::pow(rho,params->gamma));
  }

  for (int fpt=0; fpt<nFpts; fpt++) {
    auto phi = getPrimitivesFpt(fpt);
    double rho = phi[0];
    double p = phi[nDims+1];

    // Get minimum 'tau' value
    minTau = std::min(minTau, p - params->exps0*std::pow(rho,params->gamma));
  }

  if (minTau < 0) {
    // Only apply squeezing if Tau < 0; otherwise, not needed
    double rho = Uavg[0];
    double u = Uavg[1]/rho;
    double v = Uavg[2]/rho;
    double w = 0;
    if (nDims==3)
      w = Uavg[3]/rho;
    double vMagSq = u*u+v*v+w*w;
    double p = (params->gamma-1)*(Uavg[nDims+1] - 0.5*rho*vMagSq);
    double Eps = minTau / (minTau - p + params->exps0*std::pow(rho,params->gamma));

    for (int spt=0; spt<nSpts; spt++) {
      for (int i=0; i<nFields; i++) {
        U_spts(spt,i) = Eps*Uavg[i] + (1-Eps)*U_spts(spt,i);
      }
    }

    for (int fpt=0; fpt<nFpts; fpt++) {
      for (int i=0; i<nFields; i++) {
        U_fpts(fpt,i) = Eps*Uavg[i] + (1-Eps)*U_fpts(fpt,i);
      }
    }
  }
}

void ele::checkEntropyPlot()
{
  /* --- Fisrt, check if density is negative and squeeze if needed --- */
  bool negRho = false;
  double minRho = 1e15;
  double tol = 1e-10;   // Tolerance for squeezing

  for (int spt=0; spt<nSpts; spt++) {
    if (U_spts(spt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_spts(spt,0));
    }
  }

  for (int fpt=0; fpt<nFpts; fpt++) {
    if (U_fpts(fpt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_fpts(fpt,0));
    }
  }

  for (int mpt=0; mpt<nMpts; mpt++) {
    if (U_mpts(mpt,0) < 0) {
      negRho = true;
      minRho = min(minRho,U_mpts(mpt,0));
    }
  }

  // --- Do the squeezing on density (if needed) ---
  if (negRho) {
    double eps = abs(Uavg[0] - tol)/(Uavg[0] - minRho);
    for (int spt=0; spt<nSpts; spt++) {
      U_spts(spt,0) = (1-eps)*Uavg[0] + eps*U_spts(spt,0);
    }

    for (int fpt=0; fpt<nFpts; fpt++) {
      U_fpts(fpt,0) = (1-eps)*Uavg[0] + eps*U_fpts(fpt,0);
    }

    for (int mpt=0; mpt<nMpts; mpt++) {
      U_mpts(mpt,0) = (1-eps)*Uavg[0] + eps*U_mpts(mpt,0);
    }
  }

  /* --- Next, check for entropy loss and correct if needed --- */

  double minTau = 1e15; // Entropy-bounding value
  for (int spt=0; spt<nSpts; spt++) {
    auto phi = getPrimitives(spt);
    double rho = phi[0];
    double p = phi[nDims+1];

    // Get minimum 'tau' value
    minTau = std::min(minTau, p - params->exps0*std::pow(rho,params->gamma));
  }

  for (int fpt=0; fpt<nFpts; fpt++) {
    auto phi = getPrimitivesFpt(fpt);
    double rho = phi[0];
    double p = phi[nDims+1];

    // Get minimum 'tau' value
    minTau = std::min(minTau, p - params->exps0*std::pow(rho,params->gamma));
  }

  for (int mpt=0; mpt<nMpts; mpt++) {
    auto phi = getPrimitivesMpt(mpt);
    double rho = phi[0];
    double p = phi[nDims+1];

    // Get minimum 'tau' value
    minTau = std::min(minTau, p - params->exps0*std::pow(rho,params->gamma));
  }

  if (minTau < 0) {
    // Only apply squeezing if Tau < 0; otherwise, not needed
    double rho = Uavg[0];
    double u = Uavg[1]/rho;
    double v = Uavg[2]/rho;
    double w = 0;
    if (nDims==3)
      w = Uavg[3]/rho;
    double vMagSq = u*u+v*v+w*w;
    double p = (params->gamma-1)*(Uavg[nDims+1] - 0.5*rho*vMagSq);
    double Eps = minTau / (minTau - p + params->exps0*std::pow(rho,params->gamma));

    for (int spt=0; spt<nSpts; spt++) {
      for (int i=0; i<nFields; i++) {
        U_spts(spt,i) = Eps*Uavg[i] + (1-Eps)*U_spts(spt,i);
      }
    }

    for (int fpt=0; fpt<nFpts; fpt++) {
      for (int i=0; i<nFields; i++) {
        U_fpts(fpt,i) = Eps*Uavg[i] + (1-Eps)*U_fpts(fpt,i);
      }
    }

    for (int mpt=0; mpt<nMpts; mpt++) {
      for (int i=0; i<nFields; i++) {
        U_mpts(mpt,i) = Eps*Uavg[i] + (1-Eps)*U_mpts(mpt,i);
      }
    }
  }
}

vector<point> ele::getPpts(void)
{
  uint nPpts = (order+3)*(order+3);
  if (nDims == 3) nPpts *= (order+3);

  vector<point> posPpts(nPpts);
  for (uint ppt = 0; ppt < nPpts; ppt++)
    posPpts[ppt] = point(&pos_ppts(ppt,0),nDims);

    return posPpts;
}

void ele::restart(ifstream &file, input* _params, geo* _Geo)
{
  params = _params;
  Geo = _Geo;

  // Get the "<Piece _ >" line
  string str;
  getline(file,str);

  stringstream ss;
  ss.precision(15);
  string str1, str2;
  ss.str(str);
  ss >> str >> str1 >> str2;

  int nPts, nCells;

  // Find quotation marks around # of points & remove
  size_t ind = str1.find("\"");
  str1.erase(str1.begin(),str1.begin()+ind+1);
  ind = str1.find("\"");
  if (ind>10) {
    cout << "rank " << params->rank << ", ind = " << ind << endl;
    cout << "Restart-file element doesn't exist!" << endl;
    for (int spt=0; spt<nSpts; spt++)
      for (int k=0; k<nFields; k++)
        U_spts(spt,k) = 100.;
    return;
  }
  str1.erase(ind,1);

  ss.str(std::string(""));
  ss.clear();  // This is how to reset stringstreams!
  ss.str(str1);
  ss >> nPts;

  // Find quotation marks around # of cells & remove
  ind = str2.find("\"");
  str2.erase(0,ind);
  ind = str2.find("\"");
  str2.erase(ind,1);
  ss.str(""); ss.clear();
  ss.str(str2);
  ss >> nCells;  // this nCells is the number of plotting sub-cells within an actual element

  nDims = params->nDims;
  if (nDims == 2) {
    order = sqrt(nCells) - 2;
    nSpts = (order+1)*(order+1);
  }
  else if (nDims == 3) {
    order = cbrt(nCells) - 2;
    nSpts = (order+1)*(order+1)*(order+1);
  }

  matrix<double> opp_interp;
  int nSpts_final;
  if (order != params->order) {
    /* Setup inter-order interpolation operator */

    if (nDims == 2)
      nSpts_final = (params->order+1)*(params->order+1);
    else
      nSpts_final = (params->order+1)*(params->order+1)*(params->order+1);

    opp_interp.setup(nSpts_final, nSpts);

    auto loc_spts_r = getPts1D(params->sptsTypeQuad,order);
    auto loc_spts_f = getPts1D(params->sptsTypeQuad,params->order);

    point loc;
    if (nDims == 2) {
      for (uint spt = 0; spt < nSpts_final; spt++) {
        loc.x = loc_spts_f[spt%(params->order+1)];
        loc.y = loc_spts_f[floor(spt/(params->order+1))];
        for (uint rspt = 0; rspt < nSpts; rspt++) {
          uint ispt = rspt % (order+1);
          uint jspt = floor(rspt/(order+1));
          opp_interp(spt, rspt) = Lagrange(loc_spts_r,loc.x,ispt) * Lagrange(loc_spts_r,loc.y,jspt);
        }
      }
    }
    else {
      for (uint rspt = 0; rspt < nSpts; rspt++) {
        uint ksptr = rspt / ((order+1)*(order+1));
        uint jsptr = (rspt-(order+1)*(order+1)*ksptr)/(order+1);
        uint isptr = rspt - (order+1)*jsptr - (order+1)*(order+1)*ksptr;
        for (uint fspt = 0; fspt < nSpts_final; fspt++) {
          uint ksptf = fspt / ((params->order+1)*(params->order+1));
          uint jsptf = (fspt-(params->order+1)*(params->order+1)*ksptf)/(params->order+1);
          uint isptf = fspt - (params->order+1)*jsptf - (params->order+1)*(params->order+1)*ksptf;
          loc.x = loc_spts_f[isptf];
          loc.y = loc_spts_f[jsptf];
          loc.z = loc_spts_f[ksptf];

          opp_interp(fspt, rspt) = Lagrange(loc_spts_r,loc.x,isptr) * Lagrange(loc_spts_r,loc.y,jsptr) * Lagrange(loc_spts_r,loc.z,ksptr);
        }
      }
    }
  }

  // Move on to the first <DataArray>

  matrix<double>  tempU(nSpts,nDims);

  if (params->equation == NAVIER_STOKES) {
    matrix<double> tempV(nSpts,nDims);
    vector<double> tempP(nSpts);

    bool foundRho = false;
    bool foundV = false;
    bool foundP = false;

    while ( !(foundRho && foundV && foundP) ) {

      while (getline(file,str)) {
        ss.str(""); ss.clear();
        ss.str(str);
        ss >> str1;
        if (str1.compare("<DataArray")==0) {
          while (str1.find("Name=") == string::npos) {
            ss >> str1;
          }
          break;
        }
      }

      // Extract field name
      ind = str1.find("\"");
      str1.erase(0,ind+1);
      ind = str1.find("\"");
      str1.erase(ind,1);

      // Get the data line
      getline(file,str);
      ss.str(""); ss.clear(); ss.precision(15);
      ss.str(str);
      if (str1.compare("Density")==0) {
        foundRho = true;
        double tmp;
        if (nDims == 2) {
          // Skip the data at the mesh nodes and flux points along 'bottom' edge
          for (int i=0; i<2+(order+1); i++) {
            ss >> tmp;
          }

          // Get the data at the solution points
          for (int i=0; i<(order+1); i++) {
            ss >> tmp;
            for (int j=0; j<(order+1); j++) {
              ss >> U_spts(j+i*(order+1),0);
            }
            ss >> tmp;
          }
        }
        else {
          // Skip the data at the mesh nodes and flux points along 'bottom' face
          for (int i=0; i<(order+3)*(order+3); i++)
            ss >> tmp;

          // Get the data at the solution points
          for (int k=0; k<(order+1); k++) {
            // Skip front-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;

            for (int j=0; j<(order+1); j++) {
              ss >> tmp;
              for (int i=0; i<(order+1); i++) {
                ss >> U_spts(i+(order+1)*(j+(order+1)*k),0);
              }
              ss >> tmp;
            }

            // Skip back-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;
          }
        }
      }
      else if (str1.compare("Velocity")==0) {
        foundV = true;
        double tmp1, tmp2, tmp3;
        if (nDims == 2) {
          // Skip the data at the mesh nodes and flux points along 'bottom' edge
          for (int i=0; i<2+(order+1); i++) {
            ss >> tmp1 >> tmp2 >> tmp3;
          }

          // Get the data at the solution points; skip the flux points at
          // either end of each row
          for (int i=0; i<(order+1); i++) {
            ss >> tmp1 >> tmp2 >> tmp3;
            for (int j=0; j<(order+1); j++) {
              ss >> tempV(j+i*(order+1),0);
              ss >> tempV(j+i*(order+1),1);
              ss >> tmp1;
            }
            ss >> tmp1 >> tmp2 >> tmp3;
          }
        }
        else {
          // Skip the data at the mesh nodes and flux points along 'bottom' face
          for (int i=0; i<(order+3)*(order+3); i++)
            ss >> tmp1 >> tmp2 >> tmp3;

          // Get the data at the solution points
          for (int k=0; k<(order+1); k++) {
            // Skip front-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp1 >> tmp2 >> tmp3;

            for (int j=0; j<(order+1); j++) {
              ss >> tmp1 >> tmp2 >> tmp3;
              for (int i=0; i<(order+1); i++) {
                ss >> tempV(i+(order+1)*(j+(order+1)*k),0);
                ss >> tempV(i+(order+1)*(j+(order+1)*k),1);
                ss >> tempV(i+(order+1)*(j+(order+1)*k),2);
              }
              ss >> tmp1 >> tmp2 >> tmp3;
            }

            // Skip back-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp1 >> tmp2 >> tmp3;
          }
        }
      }
      else if (str1.compare("Pressure")==0) {
        foundP = true;
        double tmp;
        if (nDims == 2) {
          // Skip the data at the mesh nodes and flux points along 'bottom' edge
          for (int i=0; i<2+(order+1); i++) {
            ss >> tmp;
          }

          // Get the data at the solution points
          for (int i=0; i<(order+1); i++) {
            ss >> tmp;
            for (int j=0; j<(order+1); j++) {
              ss >> tempP[j+i*(order+1)];
            }
            ss >> tmp;
          }
        }
        else {
          // Skip the data at the mesh nodes and flux points along 'bottom' face
          for (int i=0; i<(order+3)*(order+3); i++)
            ss >> tmp;

          // Get the data at the solution points
          for (int k=0; k<(order+1); k++) {
            // Skip front-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;

            for (int j=0; j<(order+1); j++) {
              ss >> tmp;
              for (int i=0; i<(order+1); i++) {
                ss >> tempP[i+(order+1)*(j+(order+1)*k)];
              }
              ss >> tmp;
            }

            // Skip back-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;
          }
        }
      }
      else if (str1.compare("EntropyErr")==0 && params->calcEntropySensor) {
        double tmp;
        if (nDims == 2) {
          // Skip the data at the mesh nodes and flux points along 'bottom' edge
          for (int i=0; i<2+(order+1); i++) {
            ss >> tmp;
          }

          // Get the data at the solution points
          for (int i=0; i<(order+1); i++) {
            ss >> tmp;
            for (int j=0; j<(order+1); j++) {
              ss >> S_spts(j+i*(order+1));
            }
            ss >> tmp;
          }
        }
        else {
          // Skip the data at the mesh nodes and flux points along 'bottom' face
          for (int i=0; i<(order+3)*(order+3); i++)
            ss >> tmp;

          // Get the data at the solution points
          for (int k=0; k<(order+1); k++) {
            // Skip front-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;

            for (int j=0; j<(order+1); j++) {
              ss >> tmp;
              for (int i=0; i<(order+1); i++) {
                ss >> S_spts(i+(order+1)*(j+(order+1)*k));
              }
              ss >> tmp;
            }

            // Skip back-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;
          }
        }
      }
      else {
        // Not needed; ignore
      }
    }

    // Convert primitive variables to conservative variables
    for (int spt=0; spt<nSpts; spt++) {
      tempU(spt,1) = tempU(spt,0)*tempV(spt,0);
      tempU(spt,2) = tempU(spt,0)*tempV(spt,1);
      double vSq = tempV(spt,0)*tempV(spt,0)+tempV(spt,1)*tempV(spt,1);
      if (nDims == 3) {
        tempU(spt,3) = tempU(spt,0)*tempV(spt,2);
        vSq += tempV(spt,2)*tempV(spt,2);
      }
      tempU(spt,nDims+1) = tempP[spt]/(params->gamma-1) + 0.5*tempU(spt,0)*vSq;
    }
  }
  else if (params->equation == ADVECTION_DIFFUSION) {

    bool foundRho = false;
    while (!foundRho) {
      while (getline(file,str)) {
        ss.str(""); ss.clear();
        ss.str(str);
        ss >> str1;
        if (str1.compare("<DataArray")==0) {
          while (str1.find("Name=") == string::npos) {
            ss >> str1;
          }
          break;
        }
      }

      // Extract field name
      ind = str1.find("\"");
      str1.erase(0,ind+1);
      ind = str1.find("\"");
      str1.erase(ind,1);

      // Get the data line
      getline(file,str);
      ss.str(""); ss.clear(); ss.precision(15);
      ss.str(str);
      if (str1.compare("Density")==0) {
        foundRho = true;
        double tmp;
        if (nDims == 2) {
          // Skip the data at the mesh nodes and flux points along 'bottom' edge
          for (int i=0; i<2+(order+1); i++) {
            ss >> tmp;
          }

          // Get the data at the solution points
          for (int i=0; i<(order+1); i++) {
            ss >> tmp;
            for (int j=0; j<(order+1); j++) {
              ss >> tempU(j+i*(order+1),0);
            }
            ss >> tmp;
          }
        }
        else {
          // Skip the data at the mesh nodes and flux points along 'bottom' face
          for (int i=0; i<(order+3)*(order+3); i++)
            ss >> tmp;

          // Get the data at the solution points
          for (int k=0; k<(order+1); k++) {
            // Skip front-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;

            for (int j=0; j<(order+1); j++) {
              ss >> tmp;
              for (int i=0; i<(order+1); i++) {
                ss >> tempU(i+(order+1)*(j+(order+1)*k),0);
              }
              ss >> tmp;
            }

            // Skip back-face nodes / flux points
            for (int m=0; m<(order+3); m++)
              ss >> tmp;
          }
        }
      }
    }
  }

  // Move to end of this element's data
  getline(file,str);
  while(str.find("</Piece>")==string::npos)
    getline(file,str);

  if (order != params->order) {
    /* Use interpolation operator to interpolate from restart order to final order */
    matrix<double>  tempU2(nSpts_final,nFields);
    opp_interp.timesMatrix(tempU, tempU2);

    nSpts = nSpts_final;
    order = params->order;
    tempU.setup(nSpts,nFields);
    tempU = tempU2;
  }

  /* Copy U into global storage */

  for (uint spt = 0; spt < nSpts; spt++)
    for (uint k = 0; k < nFields; k++)
      U_spts(spt, k) = tempU(spt, k);
}

vector<double> ele::getNormResidual(int normType)
{
  vector<double> res(nFields,0);

  // Integrating residual over element using Gaussian integration
  auto weights = getQptWeights(order,nDims);

  for (int spt=0; spt<nSpts; spt++) {
    for (int i=0; i<nFields; i++) {
      if (normType == 1) {
        res[i] += abs(divF_spts(0,spt,i)) * weights[spt];
      }
      else if (normType == 2) {
        res[i] += divF_spts(0,spt,i)*divF_spts(0,spt,i) / detJac_spts(spt) * weights[spt];
      }
      else if (normType == 3) {
        // Infinity norm
        res[i] = max(abs(divF_spts(0,spt,i))/detJac_spts(spt),res[i]);
      }
    }
  }

  return res;
}

point ele::getPosSpt(uint spt)
{
  return point(&pos_spts(spt,0),nDims);
}

point ele::getPosFpt(uint fpt)
{
  return point(&pos_fpts(fpt,0),nDims);
}

void ele::getPosSpts(double* posSpts)
{
  for (int spt=0; spt<nSpts; spt++)
    for (int dim=0; dim<nDims; dim++)
      posSpts[spt*nDims+dim] = pos_spts(spt,dim);
}

vector<point> ele::getPosSpts(void)
{
  vector<point> posSpts(nSpts);

  for (int spt=0; spt<nSpts; spt++)
    posSpts[spt] = point(&pos_spts(spt,0),nDims);

  return posSpts;
}

uint ele::getNDims() const
{
  return nDims;
}

void ele::setNDims(int value)
{
  nDims = value;
}
uint ele::getNFields() const
{
  return nFields;
}

void ele::setNFields(int value)
{
  nFields = value;
}
uint ele::getNSpts() const
{
  return nSpts;
}

void ele::setNSpts(int value)
{
  nSpts = value;
}
uint ele::getNFpts() const
{
  return nFpts;
}

void ele::setNFpts(int value)
{
  nFpts = value;
}

double ele::getSensor(void)
{
   return sensor;
}

void ele::getUSpts(double* Uvec)
{
  for (int spt=0; spt<nSpts; spt++)
    for (int field=0; field<nFields; field++)
      Uvec[spt*nFields+field] = U_spts(spt,field);
}

void ele::setUSpts(double* Uvec)
{
  for (int spt=0; spt<nSpts; spt++)
    for (int field=0; field<nFields; field++)
      U_spts(spt,field) = Uvec[spt*nFields+field];
}

