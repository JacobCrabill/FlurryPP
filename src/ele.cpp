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
  Geo = in_Geo;

  nodes.clear();
  for (auto &n: in_nodes)
    nodes.push_back(n);

}

void ele::initialize(void)
{

}

void ele::setup(input *inParams, geo *inGeo)
{
  /* --- Basic Stuff --- */
  params = inParams;
  Geo = inGeo;

  order = params->order;
  nDims = params->nDims;
  nFields = params->nFields;

  if (eType == QUAD || eType == HEX)
    sptsType = params->sptsTypeQuad;
  else
    FatalError("Only quads and hexes implemented.");

  loc_spts = getLocSpts(eType,order,sptsType);
  loc_fpts = getLocFpts(eType,order,sptsType);

  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  pos_spts.resize(nSpts);
  pos_fpts.resize(nFpts);

  /* --- Setup all data arrays --- */
  setupArrays();

  /* --- Final Step: calculate physical->reference transforms
   * and store shape basis values for future use --- */
  setupAllGeometry();

}

void ele::setupArrays(void)
{
  U_spts.setup(nSpts,nFields);
  U_fpts.setup(nFpts,nFields);
  U_mpts.setup(nMpts,nFields);
  disFn_fpts.setup(nFpts,nFields);
  dFn_fpts.setup(nFpts,nFields);
  Fn_fpts.setup(nFpts,nFields);
  Fn_fpts.initializeToZero();

  nRKSteps = params->nRKSteps;
  divF_spts.resize(nRKSteps);
  for (auto& dF:divF_spts) dF.setup(nSpts,nFields);

  if (params->motion || params->viscous) {
    dU_spts.resize(nDims);
    dU_fpts.resize(nDims);
    for (int dim=0; dim<nDims; dim++) {
      dU_spts[dim].setup(nSpts,nFields);
      dU_fpts[dim].setup(nFpts,nFields);
      dU_spts[dim].initializeToZero();
      dU_fpts[dim].initializeToZero();
    }
  }

  F_spts.resize(nDims);
  F_fpts.resize(nDims);
  dF_spts.setup(nDims,nDims);
  tdF_spts.resize(nDims);
  for (int i=0; i<nDims; i++) {
    F_spts[i].setup(nSpts,nFields);
    F_fpts[i].setup(nFpts,nFields);
    tdF_spts[i].setup(nSpts,nFields);
    for (int j=0; j<nDims; j++) {
      dF_spts(i,j).setup(nSpts,nFields);
    }
  }

  detJac_spts.resize(nSpts);
  detJac_fpts.resize(nFpts);
  Jac_spts.resize(nSpts);
  Jac_fpts.resize(nFpts);
  JGinv_spts.resize(nSpts);
  JGinv_fpts.resize(nFpts);
  for (auto& spt:Jac_spts) spt.setup(nDims,nDims);
  for (auto& fpt:Jac_fpts) fpt.setup(nDims,nDims);
  for (auto& spt:JGinv_spts) spt.setup(nDims,nDims);
  for (auto& fpt:JGinv_fpts) fpt.setup(nDims,nDims);

  norm_fpts.setup(nFpts,nDims);
  tNorm_fpts.setup(nFpts,nDims);
  dA_fpts.resize(nFpts);
  waveSp_fpts.resize(nFpts);

  gridVel_nodes.setup(nMpts,nDims);
  gridVel_spts.setup(nSpts,nDims);
  gridVel_fpts.setup(nFpts,nDims);

  if (params->motion != 0) {
    nodesRK = nodes;
  }

  if (params->viscous) {
    Uc_fpts.setup(nFpts,nFields);
    Uc_fpts.initializeToZero();
    dUc_fpts.setup(nFpts,nFields);
    dUc_fpts.initializeToZero();
    dU_fpts.resize(nDims);
    for (auto &du:dU_fpts) du.setup(nFpts,nFields);
  }

  if (params->scFlag)
    sensor = 0;

  if (params->equation == NAVIER_STOKES && params->calcEntropySensor) {
    S_spts.setup(nSpts,1);
    S_fpts.setup(nFpts,1);
    S_mpts.setup(nMpts,1);
  }

  tempF.setup(nDims,nFields);
  tempU.assign(nFields,0);
}

void ele::setupAllGeometry(void) {
  setShape_spts();
  setShape_fpts();
  setDShape_spts();
  setDShape_fpts();
  setTransformedNormals_fpts();
  calcTransforms();

  calcPosSpts();
  calcPosFpts();
  setPpts();
}

void ele::move(bool doTransforms)
{
  for (int i=0; i<nNodes; i++) {
    nodesRK[i] = point(Geo->xv[Geo->c2v(ID,i)],nDims);
  }

  if (params->meshType == OVERSET_MESH) {
    // Only needed for overset connectivity purposes [also called from output.cpp]
    updatePosSpts();
    updatePosFpts();
  }

  if (doTransforms) {
    calcTransforms(true);
    calcGridVelocity();
  }
}

void ele::calcGridVelocity(void)
{
  for (int iv=0; iv<nNodes; iv++) {
    for (int dim=0; dim<nDims; dim++) {
      gridVel_nodes(iv,dim) = Geo->gridVel(Geo->c2v(ID,iv),dim);
    }
  }

  gridVel_spts.initializeToZero();
  for (int spt=0; spt<nSpts; spt++) {
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        gridVel_spts(spt,dim) += shape_spts(spt,iv)*gridVel_nodes(iv,dim);
      }
    }
  }

  gridVel_fpts.initializeToZero();
  for (int fpt=0; fpt<nFpts; fpt++) {
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        gridVel_fpts(fpt,dim) += shape_fpts(fpt,iv)*gridVel_nodes(iv,dim);
      }
    }
  }
}

void ele::setShape_spts(void)
{
  shape_spts.setup(nSpts,nNodes);

  for (int spt=0; spt<nSpts; spt++) {
    switch(eType) {
      case TRI:
        shape_tri(loc_spts[spt],shape_spts[spt]);
        break;
      case QUAD:
        shape_quad(loc_spts[spt],shape_spts[spt],nNodes);
        break;
      case HEX:
        shape_hex(loc_spts[spt],shape_spts[spt],nNodes);
        break;
    }
  }
}

void ele::setShape_fpts(void)
{
  shape_fpts.setup(nFpts,nNodes);

  for (int fpt=0; fpt<nFpts; fpt++) {
    switch(eType) {
      case TRI:
        shape_tri(loc_fpts[fpt],shape_fpts[fpt]);
        break;
      case QUAD:
        shape_quad(loc_fpts[fpt],shape_fpts[fpt],nNodes);
        break;
      case HEX:
        shape_hex(loc_fpts[fpt],shape_fpts[fpt],nNodes);
        break;
    }
  }
}

void ele::setDShape_spts(void)
{
  dShape_spts.resize(nSpts);
  for (auto& dS:dShape_spts) dS.setup(nNodes,nDims);

  for (int spt=0; spt<nSpts; spt++) {
    switch(eType) {
      case TRI:
        dshape_tri(loc_spts[spt], dShape_spts[spt]);
        break;
      case QUAD:
        dshape_quad(loc_spts[spt], dShape_spts[spt],nNodes);
        break;
      case HEX:
        dshape_hex(loc_spts[spt], dShape_spts[spt],nNodes);
        break;
      default:
        FatalError("Element type not yet implemented.")
    }
  }
}

void ele::setDShape_fpts(void)
{
  dShape_fpts.resize(nFpts);
  for (auto& dS:dShape_fpts) dS.setup(nNodes,nDims);

  for (int fpt=0; fpt<nFpts; fpt++) {
    switch(eType) {
      case TRI:
        dshape_tri(loc_fpts[fpt], dShape_fpts[fpt]);
        break;
      case QUAD:
        dshape_quad(loc_fpts[fpt], dShape_fpts[fpt],nNodes);
        break;
      case HEX:
        dshape_hex(loc_fpts[fpt], dShape_fpts[fpt],nNodes);
        break;
      default:
        FatalError("Element type not yet implemented.")
    }
  }
}

void ele::setTransformedNormals_fpts(void)
{
  // Setting unit normal vector in the parent domain
  for (int fpt=0; fpt<nFpts; fpt++) {
    uint iFace;
    // Calculate shape derivatives [in the future, should pre-calculate & store]
    switch(eType) {
      case TRI:
        iFace = floor(fpt / (order+1));
        switch(iFace) {
          case 0:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = -1;
            break;
          case 1:
            tNorm_fpts(fpt,0) = sqrt(2);
            tNorm_fpts(fpt,1) = sqrt(2);
            break;
          case 2:
            tNorm_fpts(fpt,0) = -1;
            tNorm_fpts(fpt,1) = 0;
            break;
        }
        break;

      case QUAD:
        iFace = floor(fpt / (order+1));
        // Face ordering for quads: Bottom, Right, Top, Left
        switch(iFace) {
          case 0:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = -1;
            break;
          case 1:
            tNorm_fpts(fpt,0) = 1;
            tNorm_fpts(fpt,1) = 0;
            break;
          case 2:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = 1;
            break;
          case 3:
            tNorm_fpts(fpt,0) = -1;
            tNorm_fpts(fpt,1) = 0;
            break;
        }
        break;

      case HEX:
        iFace = floor(fpt / ((order+1)*(order+1)));
        switch(iFace) {
          case 0:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = 0;
            tNorm_fpts(fpt,2) = -1;
            break;
          case 1:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = 0;
            tNorm_fpts(fpt,2) = 1;
            break;
          case 2:
            tNorm_fpts(fpt,0) = -1;
            tNorm_fpts(fpt,1) = 0;
            tNorm_fpts(fpt,2) = 0;
            break;
          case 3:
            tNorm_fpts(fpt,0) = 1;
            tNorm_fpts(fpt,1) = 0;
            tNorm_fpts(fpt,2) = 0;
            break;
          case 4:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = -1;
            tNorm_fpts(fpt,2) = 0;
            break;
          case 5:
            tNorm_fpts(fpt,0) = 0;
            tNorm_fpts(fpt,1) = 1;
            tNorm_fpts(fpt,2) = 0;
            break;
        }
        break;

      default:
        FatalError("Element type not yet implemented.")
    }
  }
}

void ele::calcTransforms(bool moving)
{
  /* --- Calculate Transformation at Solution Points --- */
  for (int spt=0; spt<nSpts; spt++) {
    Jac_spts[spt].initializeToZero();

    if (!moving) {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_spts[spt](dim1,dim2) += dShape_spts[spt](i,dim2)*nodes[i][dim1];
    }
    else {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_spts[spt](dim1,dim2) += dShape_spts[spt](i,dim2)*nodesRK[i][dim1];
    }


    if (nDims == 2) {
      // Determinant of transformation matrix
      detJac_spts[spt] = Jac_spts[spt](0,0)*Jac_spts[spt](1,1)-Jac_spts[spt](1,0)*Jac_spts[spt](0,1);
      // Inverse of transformation matrix (times its determinant)
      JGinv_spts[spt](0,0) = Jac_spts[spt](1,1);  JGinv_spts[spt](0,1) =-Jac_spts[spt](0,1);
      JGinv_spts[spt](1,0) =-Jac_spts[spt](1,0);  JGinv_spts[spt](1,1) = Jac_spts[spt](0,0);
    }
    else if (nDims == 3) {
      double xr = Jac_spts[spt](0,0);   double xs = Jac_spts[spt](0,1);   double xt = Jac_spts[spt](0,2);
      double yr = Jac_spts[spt](1,0);   double ys = Jac_spts[spt](1,1);   double yt = Jac_spts[spt](1,2);
      double zr = Jac_spts[spt](2,0);   double zs = Jac_spts[spt](2,1);   double zt = Jac_spts[spt](2,2);
      detJac_spts[spt] = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);

      JGinv_spts[spt](0,0) = ys*zt - yt*zs;  JGinv_spts[spt](0,1) = xt*zs - xs*zt;  JGinv_spts[spt](0,2) = xs*yt - xt*ys;
      JGinv_spts[spt](1,0) = yt*zr - yr*zt;  JGinv_spts[spt](1,1) = xr*zt - xt*zr;  JGinv_spts[spt](1,2) = xt*yr - xr*yt;
      JGinv_spts[spt](2,0) = yr*zs - ys*zr;  JGinv_spts[spt](2,1) = xs*zr - xr*zs;  JGinv_spts[spt](2,2) = xr*ys - xs*yr;
    }
    if (detJac_spts[spt]<0) FatalError("Negative Jacobian at solution points.");
  }

  /* --- Calculate Transformation at Flux Points --- */
  for (int fpt=0; fpt<nFpts; fpt++) {
    // Calculate transformation Jacobian matrix - [dx/dr, dx/ds; dy/dr, dy/ds]
    Jac_fpts[fpt].initializeToZero();
    if (!moving) {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_fpts[fpt][dim1][dim2] += dShape_fpts[fpt][i][dim2]*nodes[i][dim1];
    }
    else {
      for (int i=0; i<nNodes; i++)
        for (int dim1=0; dim1<nDims; dim1++)
          for (int dim2=0; dim2<nDims; dim2++)
            Jac_fpts[fpt][dim1][dim2] += dShape_fpts[fpt][i][dim2]*nodesRK[i][dim1];
    }


    if (nDims == 2) {
      detJac_fpts[fpt] = Jac_fpts[fpt](0,0)*Jac_fpts[fpt](1,1)-Jac_fpts[fpt](1,0)*Jac_fpts[fpt](0,1);
      // Inverse of transformation matrix (times its determinant)
      JGinv_fpts[fpt](0,0) = Jac_fpts[fpt](1,1);  JGinv_fpts[fpt](0,1) =-Jac_fpts[fpt](0,1);
      JGinv_fpts[fpt](1,0) =-Jac_fpts[fpt](1,0);  JGinv_fpts[fpt](1,1) = Jac_fpts[fpt](0,0);
    }
    else if (nDims == 3) {
      double xr = Jac_fpts[fpt](0,0);   double xs = Jac_fpts[fpt](0,1);   double xt = Jac_fpts[fpt](0,2);
      double yr = Jac_fpts[fpt](1,0);   double ys = Jac_fpts[fpt](1,1);   double yt = Jac_fpts[fpt](1,2);
      double zr = Jac_fpts[fpt](2,0);   double zs = Jac_fpts[fpt](2,1);   double zt = Jac_fpts[fpt](2,2);
      detJac_fpts[fpt] = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);
      // Inverse of transformation matrix (times its determinant)
      JGinv_fpts[fpt](0,0) = ys*zt - yt*zs;  JGinv_fpts[fpt](0,1) = xt*zs - xs*zt;  JGinv_fpts[fpt](0,2) = xs*yt - xt*ys;
      JGinv_fpts[fpt](1,0) = yt*zr - yr*zt;  JGinv_fpts[fpt](1,1) = xr*zt - xt*zr;  JGinv_fpts[fpt](1,2) = xt*yr - xr*yt;
      JGinv_fpts[fpt](2,0) = yr*zs - ys*zr;  JGinv_fpts[fpt](2,1) = xs*zr - xr*zs;  JGinv_fpts[fpt](2,2) = xr*ys - xs*yr;
    }

    /* --- Calculate outward unit normal vector at flux point --- */
    // Transform face normal from reference to physical space [JGinv .dot. tNorm]
    for (int dim1=0; dim1<nDims; dim1++) {
      norm_fpts(fpt,dim1) = 0.;
      for (int dim2=0; dim2<nDims; dim2++) {
        norm_fpts(fpt,dim1) += JGinv_fpts[fpt](dim2,dim1) * tNorm_fpts(fpt,dim2);
      }
    }

    // Store magnitude of face normal (equivalent to face area in finite-volume land)
    dA_fpts[fpt] = 0;
    for (int dim=0; dim<nDims; dim++)
      dA_fpts[fpt] += norm_fpts(fpt,dim)*norm_fpts(fpt,dim);
    dA_fpts[fpt] = sqrt(dA_fpts[fpt]);

    // Normalize
    // If we have a collapsed edge, the dA will be 0, so just set the normal to 0
    // (A normal vector at a point doesn't make sense anyways)
    if (std::fabs(dA_fpts[fpt]) < 1e-10) {
      dA_fpts[fpt] = 0.;
      for (int dim=0; dim<nDims; dim++)
        norm_fpts(fpt,dim) = 0;
    }
    else {
      for (int dim=0; dim<nDims; dim++)
        norm_fpts(fpt,dim) /= dA_fpts[fpt];
    }
  }
}

void ele::calcTransforms_point(matrix<double> &jacobian, matrix<double> &JGinv, double &detJac, const point &loc)
{
  jacobian.setup(nDims,nDims);
  jacobian.initializeToZero();
  JGinv.setup(nDims,nDims);
  JGinv.initializeToZero();

  matrix<double> dshape;
  if (nDims==2)
    dshape_quad(loc, dshape, nNodes);
  else
    dshape_hex(loc, dshape, nNodes);

  if (!params->motion) {
    for (int i=0; i<nNodes; i++)
      for (int dim1=0; dim1<nDims; dim1++)
        for (int dim2=0; dim2<nDims; dim2++)
          jacobian(dim1,dim2) += dshape(i,dim2)*nodes[i][dim1];
  }
  else {
    for (int i=0; i<nNodes; i++)
      for (int dim1=0; dim1<nDims; dim1++)
        for (int dim2=0; dim2<nDims; dim2++)
          jacobian(dim1,dim2) += dshape(i,dim2)*nodesRK[i][dim1];
  }

  if (nDims == 2) {
    // Determinant of transformation matrix
    detJac = jacobian(0,0)*jacobian(1,1)-jacobian(1,0)*jacobian(0,1);
    // Inverse of transformation matrix (times its determinant)
    JGinv(0,0) = jacobian(1,1);  JGinv(0,1) =-jacobian(0,1);
    JGinv(1,0) =-jacobian(1,0);  JGinv(1,1) = jacobian(0,0);
  } else {
    double xr = jacobian(0,0);   double xs = jacobian(0,1);   double xt = jacobian(0,2);
    double yr = jacobian(1,0);   double ys = jacobian(1,1);   double yt = jacobian(1,2);
    double zr = jacobian(2,0);   double zs = jacobian(2,1);   double zt = jacobian(2,2);
    detJac = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);

    JGinv(0,0) = ys*zt - yt*zs;  JGinv(0,1) = xt*zs - xs*zt;  JGinv(0,2) = xs*yt - xt*ys;
    JGinv(1,0) = yt*zr - yr*zt;  JGinv(1,1) = xr*zt - xt*zr;  JGinv(1,2) = xt*yr - xr*yt;
    JGinv(2,0) = yr*zs - ys*zr;  JGinv(2,1) = xs*zr - xr*zs;  JGinv(2,2) = xr*ys - xs*yr;
  }
  if (detJac<0) FatalError("Negative Jacobian at given point.");
}

point ele::calcPos(const point &loc)
{
  vector<double> shape;
  getShape(loc,shape);

  point pt;
  if (params->motion == 0) {
    for (int iv=0; iv<nNodes; iv++)
      for (int dim=0; dim<nDims; dim++)
        pt[dim] += shape[iv]*nodes[iv][dim];
  }
  else {
    for (int iv=0; iv<nNodes; iv++)
      for (int dim=0; dim<nDims; dim++)
        pt[dim] += shape[iv]*nodesRK[iv][dim];
  }

  return pt;
}

vector<double> ele::getBoundingBox(void)
{
  vector<double> bbox = {INFINITY,INFINITY,INFINITY,-INFINITY,-INFINITY,-INFINITY};
  if (params->motion == 0) {
    for (auto &pt:nodes) {
      for (int dim=0; dim<3; dim++) {
        bbox[dim]   = min(bbox[dim],  pt[dim]);
        bbox[dim+3] = max(bbox[dim+3],pt[dim]);
      }
    }
  }
  else {
    for (auto &pt:nodesRK) {
      for (int dim=0; dim<3; dim++) {
        bbox[dim]   = min(bbox[dim],  pt[dim]);
        bbox[dim+3] = max(bbox[dim+3],pt[dim]);
      }
    }
  }

  return bbox;
}

point ele::getRefLoc(const point &pos)
{
  // --- NOTE: Need a better method for high-aspect-ratio elements!! ---
  // Tested in Matlab; fails if element if too thin in one direction, even if
  // point is within element
  // What about using Nelder-Meade to minimiz norm(dx)?
  point xin;  // Current iterate xi_n (in reference space)
  point xn;   // Current iterate's physical position
  Vec3 dx;    // Difference from current iterate to desired position

  matrix<double> J(nDims,nDims);
  matrix<double> Jinv(nDims,nDims);

  double tol = 1e-6;  // Want to go to smaller tolerance once all working

  // No direct inverse-isoparametric mapping, so must iterate (Newton's Method)
  // Note: Initial guess is (0,0,0) by point()
  xn = calcPos(xin);
  dx = xn - pos;

  int iter = 0;
  while (dx.norm() > tol && iter < 200) {
    xin -= Jinv*dx;
    getInverseMapping(xin,J,Jinv);
    xn = calcPos(xin);
    dx = xn - pos;
    iter++;
  }

  return xin;
}

void ele::getInverseMapping(const point xi, matrix<double> &J, matrix<double> &Jinv)
{
  matrix<double> dshape(nNodes,nDims);
  J.setup(nDims,nDims);    J.initializeToZero();
  Jinv.setup(nDims,nDims); Jinv.initializeToZero();
  double detJ = 0;

  // Get shape (isoparametric mapping) derivatives at mesh nodes
  switch (eType) {
  case QUAD:
    dshape_quad(xi,dshape,nNodes);
    break;

  case HEX:
    dshape_hex(xi,dshape,nNodes);
    break;
  }

  // Calculate Jacobian matrix and its inverse
  if (nDims == 2) {
    for (int n=0; n<nNodes; n++)
      for (int dim2=0; dim2<nDims; dim2++)
        for (int dim1=0; dim1<nDims; dim1++)
          J(dim2,dim1) += dshape(n,dim2)*nodes[n][dim1];

    detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    if (detJ <= 0) FatalError("Negative Jacobian in inverse-mapping calculation");

    Jinv(0,0) = J(1,1)/detJ;  Jinv(0,1) =-J(0,1)/detJ;
    Jinv(1,0) =-J(1,0)/detJ;  Jinv(1,1) = J(0,0)/detJ;
  }
  else {
    for (int n=0; n<nNodes; n++)
      for (int dim2=0; dim2<nDims; dim2++)
        for (int dim1=0; dim1<nDims; dim1++)
          J(dim2,dim1) += dshape(n,dim2)*nodes[n][dim1];

    double xr = J(0,0);   double xs = J(0,1);   double xt = J(0,2);
    double yr = J(1,0);   double ys = J(1,1);   double yt = J(1,2);
    double zr = J(2,0);   double zs = J(2,1);   double zt = J(2,2);
    detJ = xr*(ys*zt - yt*zs) - xs*(yr*zt - yt*zr) + xt*(yr*zs - ys*zr);
    if (detJ <= 0) FatalError("Negative Jacobian in inverse-mapping calculation");

        Jinv(0,0) = ys*zt - yt*zs;  Jinv(0,1) = xt*zs - xs*zt;  Jinv(0,2) = xs*yt - xt*ys;
    Jinv(1,0) = yt*zr - yr*zt;  Jinv(1,1) = xr*zt - xt*zr;  Jinv(1,2) = xt*yr - xr*yt;
    Jinv(2,0) = yr*zs - ys*zr;  Jinv(2,1) = xs*zr - xr*zs;  Jinv(2,2) = xr*ys - xs*yr;

    for (int i=0; i<nDims; i++)
      for (int j=0; j<nDims; j++)
        Jinv(i,j) /= detJ;
  }
}

double ele::getDxNelderMeade(point refLoc, point physPos)
{
  point pt = calcPos(refLoc);
  Vec3 dx = physPos - pt;
  return dx.norm();
}

bool ele::getRefLocNelderMeade(point pos, point& loc)
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
  vector<double> F(nPts);
  matrix<double> X(nPts,nVars);

  // Starting location for search
  double L = .75;
  if (nDims == 3) {
    X(0,0) =-L*.5; X(0,1) =-L*.43301; X(0,2) =-L*.375;
    X(1,0) = L*.5; X(1,1) =-L*.43301; X(1,2) =-L*.375;
    X(2,0) = L*0;  X(2,1) = L*.43301; X(2,2) =-L*.375;
    X(3,0) = L*0;  X(3,1) = L*0;      X(3,2) = L*.375;
  }
  else {
    X(0,0) =-L*.5; X(0,1) =-L*.43301;
    X(1,0) = L*.5; X(1,1) =-L*.43301;
    X(2,0) = L*0;  X(2,1) = L*.43301;
  }

  // Evaluate the 'function' at the initial 'points'
  for (int i=0; i<nPts; i++)
    F[i] = getDxNelderMeade(point(X[i],nDims),pos);

  double tol = 1e-11;
  int iter = 0;
  while (iter < 300 && getMin(F)>tol) {
    auto ind = getOrder(F);
    point Xn = point(X[ind[nPts-1]],nDims);  // Point with the highest value of F
    point X0;                          // Centroid of all other points
    point Xr;                          // Reflected point

    // Take centroid of all points besides Xn
    for (int j=0; j<nPts-1; j++)
      X0 += point(X[ind[j]],nDims)/(nPts-1);
    // Reflect Xn around X0
    Xr = X0 + (X0-Xn);

    double Fr = getDxNelderMeade(Xr,pos);

    // Determine what to do with the new point
    if (Fr < F[ind[nPts-2]]) {
      // We will be keeping this point
      if (Fr < F[ind[0]]) {
        // This one's good; keep going! Expand from Xr
        point Xe = Xr + (X0-Xn);
        double Fe = getDxNelderMeade(Xe,pos);

        if (Fe < Fr) {
          // This one's even better; use it instead
          for (int i=0; i<nVars; i++) {
            X(ind[nPts-1],i) = Xe[i];
            F[ind[nPts-1]] = Fe;
          }
        }
        else {
          // Xe/Fe was no better; stick with Fr, Xr
          for (int i=0; i<nVars; i++) {
            X(ind[nPts-1],i) = Xr[i];
            F[ind[nPts-1]] = Fr;
          }
        }
      }
      else {
        // This one's somewhere in the middle; replace Xn with Xr
        for (int i=0; i<nVars; i++) {
          X(ind[nPts-1],i) = Xr[i];
          F[ind[nPts-1]] = Fr;
        }
      }
    }
    else {
      // Try reducing the size of the simplex
      point Xc = X0 - (X0-Xn)*.5;
      double Fc = getDxNelderMeade(Xc,pos);
      if (Fc < F[ind[nPts-1]]) {
        // Bringing this point in is better; use it
        for (int i=0; i<nVars; i++) {
          X(ind[nPts-1],i) = Xc[i];
          F[ind[nPts-1]] = Fc;
        }
      }
      else {
        // Bringing this point in didn't work; shrink the simplex onto
        // the smallest-valued vertex
        point X1 = point(X[ind[0]],nDims);
        for (int i=1; i<nPts; i++) {
          for (int j=0; j<nVars; j++) {
            X(ind[i],j) = X1[j] + 0.5*(X(ind[i],j)-X1[j]);
          }
          point xTmp = point(X[ind[i]],nDims);
          F[ind[i]] = getDxNelderMeade(xTmp,pos);
        }
      }
    }

    // Continue to iterate
    iter++;
  }

  auto ind = getOrder(F);
  loc = point(X[ind[nPts-1]],nDims);

  // Check to see if final location lies within element or not
  eps = 1e-8;
  if (std::abs(loc.x)-eps<=1 && std::abs(loc.y)-eps<=1 && std::abs(loc.z)-eps<=1 && !std::isnan(loc.norm()))
    return true;
  else
    return false;
}

//bool ele::getRefLocQuad(point pos, point &loc)
//{
//  // Check if pos lies within ele using ray tracing approach [LINEAR QUADS ONLY]

//  // Side 1
//  Vec3 dx1 = nodes[1] - nodes[0];
//  Vec3 dx2 = pos - nodes[0];
//  Vec3 cross = dx1.cross(dx2);
//  if (cross.z < 0)
//    return false;

//  // Side 2
//  dx1 = nodes[2] - nodes[1];
//  dx2 = pos - nodes[1];
//  cross = dx1.cross(dx2);
//  if (cross.z < 0)
//    return false;

//  // Side 3
//  dx1 = nodes[3] - nodes[2];
//  dx2 = pos - nodes[2];
//  cross = dx1.cross(dx2);
//  if (cross.z < 0)
//    return false;

//  // Side 4
//  dx1 = nodes[0] - nodes[3];
//  dx2 = pos - nodes[3];
//  cross = dx1.cross(dx2);
//  if (cross.z < 0)
//    return false;

//  // Point must lie within element
//}

void ele::calcPosSpts(void)
{
  for (int spt=0; spt<nSpts; spt++) {
    pos_spts[spt].zero();
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        pos_spts[spt][dim] += shape_spts(spt,iv)*nodes[iv][dim];
      }
    }
  }
}

void ele::calcPosFpts(void)
{
  for (int fpt=0; fpt<nFpts; fpt++) {
    pos_fpts[fpt].zero();
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        pos_fpts[fpt][dim] += shape_fpts(fpt,iv)*nodes[iv][dim];
      }
    }
  }
}

void ele::updatePosSpts(void)
{
  for (int spt=0; spt<nSpts; spt++) {
    pos_spts[spt].zero();
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        pos_spts[spt][dim] += shape_spts(spt,iv)*nodesRK[iv][dim];
      }
    }
  }
}

void ele::updatePosFpts(void)
{
  for (int fpt=0; fpt<nFpts; fpt++) {
    pos_fpts[fpt].zero();
    for (int iv=0; iv<nNodes; iv++) {
      for (int dim=0; dim<nDims; dim++) {
        pos_fpts[fpt][dim] += shape_fpts(fpt,iv)*nodesRK[iv][dim];
      }
    }
  }
}


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
        double x = pos_spts[spt].x;
        double y = pos_spts[spt].y;

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
        double x = pos_spts[spt].x;
        double y = pos_spts[spt].y;

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
        double r2 = pos_spts[spt]*pos_spts[spt];
        U_spts(spt,0) = exp(-r2);
      }
    }
    else if (params->icType == 1) {
      /* --- Test case for debugging - linear solution x+y+z over domain --- */
      for (int spt=0; spt<nSpts; spt++) {
        U_spts(spt,0) = pos_spts[spt].x + pos_spts[spt].y + pos_spts[spt].z;
      }
    }
    else if (params->icType == 2) {
      /* --- Test case for debugging - cos(x)*cos(y)*cos(z) over domain --- */
      for (int spt=0; spt<nSpts; spt++)
        U_spts(spt,0) = cos(2*pi*pos_spts[spt].x/6.)*cos(2*pi*pos_spts[spt].y/6.)*cos(2*pi*pos_spts[spt].z/6.);
    }
  }
}

matrix<double> ele::calcError(void)
{
  matrix<double> err(nSpts,nFields);

  if (!params->testCase) return err;

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

      double xoff = fmod( (params->time - xmin), (xmax - xmin) ) + xmin;
      double yoff = fmod( (params->time - ymin), (ymax - ymin) ) + ymin;

      for (int spt=0; spt<nSpts; spt++) {
        double x = pos_spts[spt].x - xoff;
        double y = pos_spts[spt].y - yoff;

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

      double xoff = fmod( (Uinf*cos(theta)*params->time - xmin), (xmax - xmin) ) + xmin;
      double yoff = fmod( (Uinf*sin(theta)*params->time - ymin), (ymax - ymin) ) + ymin;

      for (int spt=0; spt<nSpts; spt++) {
        double x = pos_spts[spt].x - xoff;
        double y = pos_spts[spt].y - yoff;

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

    double xoff = fmod( (params->time - xmin), (xmax - xmin) ) + xmin;
    double yoff = fmod( (params->time - ymin), (ymax - ymin) ) + ymin;
    //point off(xoff,yoff,0.);

    if (params->icType == 0) {
      /* --- Simple Gaussian bump centered at (0,0) --- */
      for (int spt=0; spt<nSpts; spt++) {
        double r2 = (pos_spts[spt].x-xoff)*(pos_spts[spt].x-xoff) + (pos_spts[spt].y-yoff)*(pos_spts[spt].y-yoff);
        err(spt,0) = exp(-r2);
      }
    }
    else if (params->icType == 2) {
      /* --- Test case for debugging - cos(x)*cos(y)*cos(z) over domain --- */
      for (int spt=0; spt<nSpts; spt++)
        err(spt,0) = cos(2*pi*pos_spts[spt].x/6.)*cos(2*pi*pos_spts[spt].y/6.)*cos(2*pi*pos_spts[spt].z/6.);
    }
  }

  // NOTE: In order to do actual L1, L2, etc. error integral, supposed to
  // take the |x| or |x|^2 FIRST, before integration!!!  Oops....
  for (int spt=0; spt<nSpts; spt++)
    for (int j=0; j<nFields; j++)
      err(spt,j) = abs(U_spts(spt,j) - err(spt,j));

  if (params->errorNorm == 2)
    for (auto &val:err.data) val *= val;  // L2 norm instead of L1

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

void ele::calcInviscidFlux_spts()
{
  for (int spt=0; spt<nSpts; spt++) {

    inviscidFlux(U_spts[spt], tempF, params);

    if (params->motion) {
      /* --- Don't transform yet; that will be handled later --- */
      for (int i=0; i<nDims; i++) {
        for (int k=0; k<nFields; k++) {
          F_spts[i](spt,k) = tempF(i,k);
        }
      }
    }
    else {
      /* --- Transform back to reference domain --- */
      for (int i=0; i<nDims; i++) {
        for (int k=0; k<nFields; k++) {
          F_spts[i](spt,k) = 0.;
          for (int j=0; j<nDims; j++) {
            F_spts[i](spt,k) += JGinv_spts[spt](i,j)*tempF(j,k);
          }
        }
      }
    }

  }
}

void ele::calcViscousFlux_spts()
{
  for (int spt=0; spt<nSpts; spt++) {

    // TEMP HACK (inefficient, but will work fine)
    matrix<double> tempDU(nDims,nFields);
    for (int dim=0; dim<nDims; dim++) {
      for (int k=0; k<nFields; k++) {
        tempDU(dim,k) = dU_spts[dim](spt,k);
      }
    }

    if (params->equation == NAVIER_STOKES)
      viscousFlux(U_spts[spt], tempDU, tempF, params);
    else if (params->equation == ADVECTION_DIFFUSION)
      viscousFluxAD(tempDU, tempF, params);

    if (params->motion) {
      /* --- Don't transform yet; that will be handled later --- */
      for (int i=0; i<nDims; i++) {
        for (int k=0; k<nFields; k++) {
          F_spts[i](spt,k) += tempF(i,k);
        }
      }
    }
    else {
      /* --- Transform back to reference domain --- */
      for (int k=0; k<nFields; k++) {
        for (int i=0; i<nDims; i++) {
          for (int j=0; j<nDims; j++) {
            F_spts[i](spt,k) += JGinv_spts[spt](i,j)*tempF(j,k);
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
    if (nDims == 2) {
      for (int spt=0; spt<nSpts; spt++) {
//        double A = gridVel_spts(spt,1)*Jac_spts[spt](0,1) - gridVel_spts(spt,0)*Jac_spts[spt](1,1);
//        double B = gridVel_spts(spt,0)*Jac_spts[spt](1,0) - gridVel_spts(spt,1)*Jac_spts[spt](0,0);
        for (int k=0; k<nFields; k++) {
          outF[0](spt,k) =  F_spts[0](spt,k)*Jac_spts[spt](1,1) - F_spts[1](spt,k)*Jac_spts[spt](0,1);// + U_spts(spt,k)*A;
          outF[1](spt,k) = -F_spts[0](spt,k)*Jac_spts[spt](1,0) + F_spts[1](spt,k)*Jac_spts[spt](0,0);// + U_spts(spt,k)*B;
        }
      }
    } else {
      for (int spt=0; spt<nSpts; spt++) {
        // Build the full 4D (space+time) Jacobian matrix & its adjoint
        matrix<double> Jacobian(4,4);
        Jacobian(3,3) = 1;
        for (int i=0; i<3; i++) {
          for (int j=0; j<3; j++)
            Jacobian(i,j) = Jac_spts[spt](i,j);
          Jacobian(i,3) = 0.;//gridVel_spts(spt,i);
        }
        matrix<double> S = Jacobian.adjoint();

        for (int dim1=0; dim1<3; dim1++)
          for (int dim2=0; dim2<3; dim2++)
            for (int k=0; k<nFields; k++)
              outF[dim1](spt,k) += F_spts[dim2](spt,k)*S(dim2,dim1);

        for (int dim=0; dim<3; dim++)
          for (int k=0; k<nFields; k++)
            outF[dim](spt,k) += U_spts(spt,k)*S(dim,3);
      }
    }
  }

  return outF;
}

vector<matrix<double>> ele::transformFlux_refToPhys(void)
{
  vector<matrix<double>> outF(nDims);

  return outF;
}

void ele::transformGradF_spts(int step)
{
  // NOTE: The 1st dim of dF is the derivative, and the 2nd is the flux direction

  if (nDims == 2) {
    for (int spt=0; spt<nSpts; spt++) {
      double A = gridVel_spts(spt,1)*Jac_spts[spt](0,1) - gridVel_spts(spt,0)*Jac_spts[spt](1,1);
      double B = gridVel_spts(spt,0)*Jac_spts[spt](1,0) - gridVel_spts(spt,1)*Jac_spts[spt](0,0);
      for (int k=0; k<nFields; k++) {
        dF_spts(0,0)(spt,k) =  dF_spts(0,0)(spt,k)*Jac_spts[spt](1,1) - dF_spts(0,1)(spt,k)*Jac_spts[spt](0,1) + dU_spts[0](spt,k)*A;
        dF_spts(1,1)(spt,k) = -dF_spts(1,0)(spt,k)*Jac_spts[spt](1,0) + dF_spts(1,1)(spt,k)*Jac_spts[spt](0,0) + dU_spts[1](spt,k)*B;
        divF_spts[step](spt,k) = dF_spts(0,0)(spt,k)+dF_spts(1,1)(spt,k);
      }
    }
  } else {
    divF_spts[step].initializeToZero();
    for (int spt=0; spt<nSpts; spt++) {
      // Build the full 4D (space+time) Jacobian matrix & its adjoint
      matrix<double> Jacobian(4,4);
      Jacobian(3,3) = 1;
      for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++)
          Jacobian(i,j) = Jac_spts[spt](i,j);
        Jacobian(i,3) = gridVel_spts(spt,i);
      }
      matrix<double> S = Jacobian.adjoint();

      for (int dim1=0; dim1<3; dim1++)
        for (int dim2=0; dim2<3; dim2++)
          for (int k=0; k<nFields; k++)
            divF_spts[step](spt,k) += dF_spts(dim2,dim1)(spt,k)*S(dim2,dim1);

      for (int dim=0; dim<3; dim++)
        for (int k=0; k<nFields; k++)
          divF_spts[step](spt,k) += dU_spts[dim](spt,k)*S(dim,3);
    }
  }
}

void ele::calcDeltaFn(void)
{
  for (int fpt=0; fpt<nFpts; fpt++) {
    for (int k=0; k<nFields; k++) {
      dFn_fpts(fpt,k) = Fn_fpts(fpt,k) - disFn_fpts(fpt,k);
    }
  }
}


void ele::calcDeltaUc(void)
{
  for (int fpt=0; fpt<nFpts; fpt++) {
    for (int k=0; k<nFields; k++) {
      dUc_fpts(fpt,k) = Uc_fpts(fpt,k) - U_fpts(fpt,k);
    }
  }
}

void ele::calcEntropyErr_spts(void)
{
  for (int spt=0; spt<nSpts; spt++) {
    auto v = getEntropyVars(spt);
    S_spts(spt) = 0;
    for (int k=0; k<nFields; k++) {
      S_spts(spt) += v[k]*divF_spts[0](spt,k);
    }
    S_spts(spt) /= detJac_spts[spt];
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
      double u = params->advectVx - gridVel_fpts(fpt,0);
      double v = params->advectVy - gridVel_fpts(fpt,1);
      double w = params->advectVz - gridVel_fpts(fpt,2);
      double csq = u*u + v*v;
      if (nDims == 3)
        csq += w*w;
      waveSp_fpts[fpt] = sqrt(csq) / dA_fpts[fpt];
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
      waveSp_fpts[fpt] = (std::abs(vN - vgN) + std::sqrt(csq)) / dA_fpts[fpt];
    }
  }
}

void ele::timeStepA(int step, double rkVal)
{
  for (int spt=0; spt<nSpts; spt++) {
    for (int i=0; i<nFields; i++) {
      U_spts(spt,i) = U0(spt,i) - rkVal * params->dt*divF_spts[step](spt,i)/detJac_spts[spt];
    }
  }
}

void ele::timeStepB(int step, double rkVal)
{
  for (int spt=0; spt<nSpts; spt++) {
    for (int i=0; i<nFields; i++) {
      U_spts(spt,i) -= rkVal * params->dt*divF_spts[step](spt,i)/detJac_spts[spt];
    }
  }
}

double ele::calcDt(void)
{
  double waveSp = 0;
  for (int fpt=0; fpt<nFpts; fpt++)
    if (dA_fpts[fpt] > 0) // ignore collapsed edges
      waveSp = max(waveSp,waveSp_fpts[fpt]);

  double dt = (params->CFL) * getCFLLimit(order) * (2.0 / (waveSp+1.e-10));
  return dt;
}


void ele::copyUspts_U0(void)
{
  U0 = U_spts;
}

void ele::copyU0_Uspts(void)
{
  U_spts = U0;
}

vector<double> ele::getPrimitives(uint spt)
{
  vector<double> V(nFields);

  if (params->equation == ADVECTION_DIFFUSION) {
    V[0] = U_spts[spt][0];
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
    V[0] = U_fpts[fpt][0];
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
    V[0] = U_mpts[mpt][0];
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
  if (eType == QUAD) {
    V.setup(nSpts+nFpts+nMpts,nFields);

    // Get solution at corner points
    for (int k=0; k<nFields; k++) {
      V(0,k)                     = U_mpts(0,k);
      V(order+2,k)               = U_mpts(1,k);
      V((order+3)*(order+3)-1,k) = U_mpts(2,k);
      V((order+3)*(order+2),k)   = U_mpts(3,k);
    }

    // Get solution at flux points
    for (int i=0; i<order+1; i++) {
      for (int k=0; k<nFields; k++) {
        V(i+1,k)                     = U_fpts(i,k);               // Bottom
        V((i+1)*(order+3),k)         = U_fpts(nFpts-i-1,k);       // Left
        V((i+2)*(order+3)-1,k)       = U_fpts(order+1+i,k);       // Right
        V((order+3)*(order+2)+i+1,k) = U_fpts(3*(order+1)-i-1,k); // Top
      }
    }

    // Get solution at solution points
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        int id = (i+1)*(order+3)+j+1;
        for (int k=0; k<nFields; k++) {
          V(id,k) = U_spts(j+i*(order+1),k);
        }
      }
    }
  }
  else if (eType == HEX) {
    int nPts1D = order+3;
    int P22 = nPts1D*nPts1D;
    int nv = 8, ne = 12;

    V.setup(nPts1D*nPts1D*nPts1D,nFields);

    for (int f=0; f<nFields; f++) {
      // Get solution at corner points
      V(0,f)                = U_mpts(0,f);
      V(order+2,f)          = U_mpts(1,f);
      V(P22-1,f)            = U_mpts(2,f);
      V(nPts1D*(order+2),f) = U_mpts(3,f);

      int base = (order+2)*P22;
      V(base,f)                    = U_mpts(4,f);
      V(base + order+2,f)          = U_mpts(5,f);
      V(base + P22-1,f)            = U_mpts(6,f);
      V(base + nPts1D*(order+2),f) = U_mpts(7,f);

      // Get solution at edge points
      for (int i=0; i<order+1; i++) {
        /* --- Bottom Edges --- */
        // edge 0-1
        V(i+1,f) = U_mpts(nv+i*ne+0,f);

        // edge 0-3
        V(nPts1D*(i+1),f) = U_mpts(nv+(order-i)*ne+3,f);

        // edge 1-2
        V(nPts1D*(i+2)-1,f) = U_mpts(nv+i*ne+1,f);

        // edge 3-2
        V(nPts1D*(order+2)+i+1,f) = U_mpts(nv+(order-i)*ne+2,f);

        /* --- Top Edges --- */
        base = P22*(order+2);
        // edge 4-5
        V(base + i+1,f) = U_mpts(nv+i*ne+4,f);

        // edge 4-7
        V(base + nPts1D*(i+1),f) = U_mpts(nv+(order-i)*ne+7,f);

        // edge 5-6
        V(base + nPts1D*(i+2)-1,f) = U_mpts(nv+i*ne+5,f);

        // edge 7-6
        V(base + nPts1D*(order+2)+i+1,f) = U_mpts(nv+(order-i)*ne+6,f);

        /* --- Mid [Vertical] Egdes --- */
        base = (i+1)*P22;
        // edge 0-4
        V(base,f) = U_mpts(nv+i*ne+8,f);

        // edge 1-5
        V(base+(order+2),f) = U_mpts(nv+i*ne+9,f);

        int base2 = nPts1D*(order+2);
        // edge 3-7
        V(base+base2,f) = U_mpts(nv+i*ne+11,f);

        // edge 2-6
        V(base+base2+order+2,f) = U_mpts(nv+i*ne+10,f);
      }

      // Get solution at flux points
      int P12 = (order+1)*(order+1);
      for (int i=0; i<order+1; i++) {
        for (int j=0; j<order+1; j++) {
          int ind1 = i + j*(order+1);
          int ind2 = order - i + (order+1)*j;
          // Bottom Face
          V(nPts1D*(j+1)+i+1,f) = U_fpts(ind1,f);

          // Top Face
          V(P22*(order+2)+(j+1)*nPts1D+i+1,f) = U_fpts(P12+ind2,f);

          // Left Face
          V(P22*(j+1)+nPts1D*(i+1),f) = U_fpts(2*P12+ind1,f);

          // Right Face
          V(P22*(j+1)+nPts1D*(i+1)+order+2,f) = U_fpts(3*P12+ind2,f);

          // Front Face
          V(P22*(j+1)+i+1,f) = U_fpts(4*P12+ind2,f);

          // Back Face
          V(P22*(j+2)+i+1-nPts1D,f) = U_fpts(5*P12+ind1,f);
        }
      }

      // Get solution at solution points
      for (int k=0; k<order+1; k++) {
        for (int j=0; j<order+1; j++) {
          for (int i=0; i<order+1; i++) {
            V(i+1+nPts1D*(j+1)+(k+1)*P22,f) = U_spts(i+(order+1)*(j+(order+1)*k),f);
          }
        }
      }
    }
  }

  if (params->equation == NAVIER_STOKES) {
    // Overwriting V, so be careful of order!
    for (uint i=0; i<V.getDim0(); i++) {
      double u = V(i,1)/V(i,0);
      double v = V(i,2)/V(i,0);
      double w = 0.;
      if (nDims == 3) w = V(i,3)/V(i,0);
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
      GV(0,dim)                = gridVel_mpts(0,dim);
      GV(order+2,dim)          = gridVel_mpts(1,dim);
      GV(P22-1,dim)            = gridVel_mpts(2,dim);
      GV(nPts1D*(order+2),dim) = gridVel_mpts(3,dim);

      int base = (order+2)*P22;
      GV(base,dim)                    = gridVel_mpts(4,dim);
      GV(base + order+2,dim)          = gridVel_mpts(5,dim);
      GV(base + P22-1,dim)            = gridVel_mpts(6,dim);
      GV(base + nPts1D*(order+2),dim) = gridVel_mpts(7,dim);

      // Get solution at edge points
      for (int i=0; i<order+1; i++) {
        /* --- Bottom Edges --- */
        // edge 0-1
        GV(i+1,dim) = gridVel_mpts(nv+i*ne+0,dim);

        // edge 0-3
        GV(nPts1D*(i+1),dim) = gridVel_mpts(nv+(order-i)*ne+3,dim);

        // edge 1-2
        GV(nPts1D*(i+2)-1,dim) = gridVel_mpts(nv+i*ne+1,dim);

        // edge 3-2
        GV(nPts1D*(order+2)+i+1,dim) = gridVel_mpts(nv+(order-i)*ne+2,dim);

        /* --- Top Edges --- */
        base = P22*(order+2);
        // edge 4-5
        GV(base + i+1,dim) = gridVel_mpts(nv+i*ne+4,dim);

        // edge 4-7
        GV(base + nPts1D*(i+1),dim) = gridVel_mpts(nv+(order-i)*ne+7,dim);

        // edge 5-6
        GV(base + nPts1D*(i+2)-1,dim) = gridVel_mpts(nv+i*ne+5,dim);

        // edge 7-6
        GV(base + nPts1D*(order+2)+i+1,dim) = gridVel_mpts(nv+(order-i)*ne+6,dim);

        /* --- Mid [Vertical] Egdes --- */
        base = (i+1)*P22;
        // edge 0-4
        GV(base,dim) = gridVel_mpts(nv+i*ne+8,dim);

        // edge 1-5
        GV(base+(order+2),dim) = gridVel_mpts(nv+i*ne+9,dim);

        int base2 = nPts1D*(order+2);
        // edge 3-7
        GV(base+base2,dim) = gridVel_mpts(nv+i*ne+11,dim);

        // edge 2-6
        GV(base+base2+order+2,dim) = gridVel_mpts(nv+i*ne+10,dim);
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
  return pos_ppts;
}

void ele::setPpts(void)
{
  if (eType == QUAD) {
    int nPts1D = order+3;
    pos_ppts.resize(nPts1D*nPts1D);

    // Get mesh (corner) points
    if (params->motion != 0) {
      pos_ppts[0*nPts1D+0]               = nodesRK[0];
      pos_ppts[0*nPts1D+order+2]         = nodesRK[1];
      pos_ppts[(order+2)*nPts1D+0]       = nodesRK[3];
      pos_ppts[(order+2)*nPts1D+order+2] = nodesRK[2];
    }
    else {
      pos_ppts[0*nPts1D+0]               = nodes[0];
      pos_ppts[0*nPts1D+order+2]         = nodes[1];
      pos_ppts[(order+2)*nPts1D+0]       = nodes[3];
      pos_ppts[(order+2)*nPts1D+order+2] = nodes[2];
    }

    // Get flux points
    for (int i=0; i<order+1; i++) {
      pos_ppts[0*nPts1D+i+1]         = pos_fpts[i];                // Bottom
      pos_ppts[(i+1)*nPts1D+0]       = pos_fpts[nFpts-i-1];        // Left
      pos_ppts[(i+1)*nPts1D+order+2] = pos_fpts[order+1+i];        // Right
      pos_ppts[(order+2)*nPts1D+i+1] = pos_fpts[3*(order+1)-i-1];  // Top
    }

    // Get solution at solution points
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        pos_ppts[(i+1)*nPts1D+j+1] = pos_spts[j+i*(order+1)];
      }
    }
  }
  else if (eType == HEX) {
    int nPts1D = order+3;
    pos_ppts.resize(nPts1D*nPts1D*nPts1D);

    int P12 = (order+1)*(order+1);
    int P22 = nPts1D*nPts1D;

    // Get mesh (corner) points
    if (params->motion != 0) {
      pos_ppts[0*nPts1D+0]               = nodesRK[0];
      pos_ppts[0*nPts1D+order+2]         = nodesRK[1];
      pos_ppts[(order+2)*nPts1D+0]       = nodesRK[3];
      pos_ppts[(order+2)*nPts1D+order+2] = nodesRK[2];

      pos_ppts[(order+2)*P22 + 0*nPts1D+0]               = nodesRK[4];
      pos_ppts[(order+2)*P22 + 0*nPts1D+order+2]         = nodesRK[5];
      pos_ppts[(order+2)*P22 + (order+2)*nPts1D+0]       = nodesRK[7];
      pos_ppts[(order+2)*P22 + (order+2)*nPts1D+order+2] = nodesRK[6];
    }
    else {
      pos_ppts[0*nPts1D+0]               = nodes[0];
      pos_ppts[0*nPts1D+order+2]         = nodes[1];
      pos_ppts[(order+2)*nPts1D+0]       = nodes[3];
      pos_ppts[(order+2)*nPts1D+order+2] = nodes[2];

      pos_ppts[(order+2)*P22 + 0*nPts1D+0]               = nodes[4];
      pos_ppts[(order+2)*P22 + 0*nPts1D+order+2]         = nodes[5];
      pos_ppts[(order+2)*P22 + (order+2)*nPts1D+0]       = nodes[7];
      pos_ppts[(order+2)*P22 + (order+2)*nPts1D+order+2] = nodes[6];
    }

    // Get edge points
    auto locPts1D = getPts1D(params->sptsTypeQuad,order);
    for (int i=0; i<order+1; i++) {
      double x1 = locPts1D[i];
      point pt, loc;
      /* --- Bottom Edges --- */
      loc.x = x1;  loc.y = -1;  loc.z = -1;  // edge 0-1
      pt = calcPos(loc);
      pos_ppts[i+1] = pt;

      loc.x = -1;  loc.y = x1;  loc.z = -1;  // edge 0-3
      pt = calcPos(loc);
      pos_ppts[nPts1D*(i+1)] = pt;

      loc.x =  1;  loc.y = x1;  loc.z = -1;  // edge 1-2
      pt = calcPos(loc);
      pos_ppts[nPts1D*(i+2)-1] = pt;

      loc.x = x1;  loc.y =  1;  loc.z = -1;  // edge 3-2
      pt = calcPos(loc);
      pos_ppts[nPts1D*(order+2)+i+1] = pt;

      /* --- Top Edges --- */
      int base = P22*(order+2);
      loc.x = x1;  loc.y = -1;  loc.z =  1;  // edge 4-5
      pt = calcPos(loc);
      pos_ppts[base + i+1] = pt;

      loc.x = -1;  loc.y = x1;  loc.z =  1;  // edge 4-7
      pt = calcPos(loc);
      pos_ppts[base + nPts1D*(i+1)] = pt;

      loc.x =  1;  loc.y = x1;  loc.z =  1;  // edge 5-6
      pt = calcPos(loc);
      pos_ppts[base + nPts1D*(i+2)-1] = pt;

      loc.x = x1;  loc.y =  1;  loc.z =  1;  // edge 7-6
      pt = calcPos(loc);
      pos_ppts[base + nPts1D*(order+2)+i+1] = pt;

      /* --- Mid [Vertical] Egdes --- */
      base = (i+1)*P22;
      loc.x = -1;  loc.y = -1;  loc.z = x1;  // edge 0-4
      pt = calcPos(loc);
      pos_ppts[base] = pt;

      loc.x =  1;  loc.y = -1;  loc.z = x1;  // edge 1-5
      pt = calcPos(loc);
      pos_ppts[base+(order+2)] = pt;

      int base2 = nPts1D*(order+2);
      loc.x = -1;  loc.y =  1;  loc.z = x1;  // edge 3-7
      pt = calcPos(loc);
      pos_ppts[base+base2] = pt;

      loc.x =  1;  loc.y =  1;  loc.z = x1;  // edge 2-6
      pt = calcPos(loc);
      pos_ppts[base+base2+order+2] = pt;
    }

    // Get flux points
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        int ind1 = i + j*(order+1);
        int ind2 = order - i + (order+1)*j;
        // Bottom Face
        pos_ppts[nPts1D*(j+1)+i+1] = pos_fpts[ind1];

        // Top Face
        pos_ppts[P22*(order+2)+(j+1)*nPts1D+i+1] = pos_fpts[P12+ind2];

        // Left Face
        pos_ppts[P22*(j+1)+nPts1D*(i+1)] = pos_fpts[2*P12+ind1];

        // Right Face
        pos_ppts[P22*(j+1)+nPts1D*(i+1)+order+2] = pos_fpts[3*P12+ind2];

        // Front Face
        pos_ppts[P22*(j+1)+i+1] = pos_fpts[4*P12+ind2];

        // Back Face
        pos_ppts[P22*(j+2)+i+1-nPts1D] = pos_fpts[5*P12+ind1];
      }
    }

    // Get solution at solution points
    for (int k=0; k<order+1; k++) {
      for (int j=0; j<order+1; j++) {
        for (int i=0; i<order+1; i++) {
          pos_ppts[i+1+nPts1D*(j+1)+(k+1)*P22] = pos_spts[i+(order+1)*(j+(order+1)*k)];
        }
      }
    }
  }
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
    nFpts = 4*(order+1);
  }
  else if (nDims == 3) {
    order = cbrt(nCells) - 2;
    nSpts = (order+1)*(order+1)*(order+1);
    nFpts = 6*(order+1)*(order+1);
  }

  if (order != params->order) {
    cout << "ele order = " << order << ", input order = " << params->order << endl;
    FatalError("Cannot restart a simulation using a different polynomial order.");
  }

  if (eType == QUAD || eType == HEX)
    sptsType = params->sptsTypeQuad;
  else
    FatalError("Only quads and hexes implemented.");

  // Move on to the first <DataArray>

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
      U_spts(spt,1) = U_spts(spt,0)*tempV(spt,0);
      U_spts(spt,2) = U_spts(spt,0)*tempV(spt,1);
      double vSq = tempV(spt,0)*tempV(spt,0)+tempV(spt,1)*tempV(spt,1);
      if (nDims == 3) {
        U_spts(spt,3) = U_spts(spt,0)*tempV(spt,2);
        vSq += tempV(spt,2)*tempV(spt,2);
      }
      U_spts(spt,nDims+1) = tempP[spt]/(params->gamma-1) + 0.5*U_spts(spt,0)*vSq;
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
    }
  }

  // Move to end of this element's data
  getline(file,str);
  while(str.find("</Piece>")==string::npos)
    getline(file,str);
}

vector<double> ele::getNormResidual(int normType)
{
  vector<double> res(nFields,0);

  // Integrating residual over element using Gaussian integration
  auto weights = getQptWeights(order,nDims);

  for (int spt=0; spt<nSpts; spt++) {
    for (int i=0; i<nFields; i++) {
      if (normType == 1) {
        res[i] += abs(divF_spts[0](spt,i)) * weights[spt];
      }
      else if (normType == 2) {
        res[i] += divF_spts[0](spt,i)*divF_spts[0](spt,i) / detJac_spts[spt] * weights[spt];
      }
      else if (normType == 3) {
        // Infinity norm
        res[i] = max(abs(divF_spts[0](spt,i))/detJac_spts[spt],res[i]);
      }
    }
  }

  return res;
}

point ele::getPosSpt(uint spt)
{
  return pos_spts[spt];
}

point ele::getPosFpt(uint fpt)
{
  return pos_fpts[fpt];
}

void ele::getPosSpts(double* posSpts)
{
  for (int spt=0; spt<nSpts; spt++)
    for (int dim=0; dim<nDims; dim++)
      posSpts[spt*3+dim] = pos_spts[spt][dim];
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

