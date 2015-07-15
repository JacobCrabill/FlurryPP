/*!
 * \file operators.cpp
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Fux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/operators.hpp"

#include <cmath>

#include "../include/polynomials.hpp"

//! Binary helper operation [for use with STL algorithms]
static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}
 
void oper::setupOperators(uint eType, uint order, geo *inGeo, input *inParams)
{
  // Get access to basic data
  Geo = inGeo;
  params = inParams;

  nDims = Geo->nDims;
  nFields = params->nFields;

  this->eType = eType;
  this->order = order;

  vector<point> loc_spts = Geo->getLocSpts(eType,order);
  vector<point> loc_fpts = Geo->getLocFpts(eType,order);

  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  // Set up each operator
  setupExtrapolateSptsFpts(loc_fpts);

  setupExtrapolateSptsMpts(loc_spts);

  setupGradSpts(loc_spts);

  setupCorrection(loc_spts,loc_fpts);

  if (params->viscous) {
    setupCorrectGradU();
  }

  // Operators needed for Shock capturing
  if (params->scFlag) {
    setupVandermonde(loc_spts);

    setupSensingMatrix();

    setupFilterMatrix();
  }
}

void oper::setupExtrapolateSptsFpts(vector<point> &loc_fpts)
{
  uint spt, fpt, ispt, jspt, kspt;
  opp_spts_to_fpts.setup(nFpts,nSpts);

  for (fpt=0; fpt<nFpts; fpt++) {
    for (spt=0; spt<nSpts; spt++) {
      switch(eType) {
        case TRI:
          opp_spts_to_fpts(fpt,spt) = eval_dubiner_basis_2d(loc_fpts[fpt],spt,order);
          break;
        case QUAD: {
          vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
          // First, get the i an j ID of the spt
          ispt = spt%(nSpts/(order+1));
          jspt = floor(spt/(order+1));
          opp_spts_to_fpts(fpt,spt) = Lagrange(locSpts1D,loc_fpts[fpt].x,ispt) * Lagrange(locSpts1D,loc_fpts[fpt].y,jspt);
          break;
        }
        case HEX: {
          vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
          // First, get the i an j ID of the spt
          kspt = spt/((order+1)*(order+1));
          jspt = (spt-(order+1)*(order+1)*kspt)/(order+1);
          ispt = spt - (order+1)*jspt - (order+1)*(order+1)*kspt;
          opp_spts_to_fpts(fpt,spt) = Lagrange(locSpts1D,loc_fpts[fpt].x,ispt) * Lagrange(locSpts1D,loc_fpts[fpt].y,jspt) * Lagrange(locSpts1D,loc_fpts[fpt].z,kspt);
          break;
        }
        default:
          FatalError("Element type not yet supported.");
      }
    }
  }
}

void oper::setupExtrapolateSptsMpts(vector<point> &loc_spts)
{
  uint nSpts = loc_spts.size();

  switch(eType) {
    case TRI: {
      opp_spts_to_mpts.setup(3,nSpts);
      point vert1, vert2, vert3;
      vert1.x = -1;  vert1.y = -1;
      vert2.x =  1;  vert2.y = -1;
      vert3.x = -1;  vert3.y =  1;
      for (uint spt=0; spt<nSpts; spt++) {
        opp_spts_to_mpts(0,spt) = eval_dubiner_basis_2d(vert1,spt,order);
        opp_spts_to_mpts(1,spt) = eval_dubiner_basis_2d(vert2,spt,order);
        opp_spts_to_mpts(2,spt) = eval_dubiner_basis_2d(vert3,spt,order);
      }
      break;
    }
    case QUAD: {
      uint ispt, jspt;
      opp_spts_to_mpts.setup(4,nSpts);
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint spt=0; spt<nSpts; spt++) {
        // First, get the i an j ID of the spt
        ispt = spt%(nSpts/(order+1));
        jspt = floor(spt/(order+1));
        // Next, get evaluate Lagrange solution basis at corners
        opp_spts_to_mpts(0,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,-1,jspt);
        opp_spts_to_mpts(1,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,-1,jspt);
        opp_spts_to_mpts(2,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D, 1,jspt);
        opp_spts_to_mpts(3,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D, 1,jspt);
      }
      break;
    }
    case HEX: {
      uint ispt, jspt, kspt;
      int nv = 8, ne = 12;
      // Have to put extra points on the edges to connect fpts & spts to mpts
      opp_spts_to_mpts.setup(nv+ne*(order+1),nSpts);
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint spt=0; spt<nSpts; spt++) {
        // First, get the i an j ID of the spt
        kspt = spt/((order+1)*(order+1));
        jspt = (spt-(order+1)*(order+1)*kspt)/(order+1);
        ispt = spt - (order+1)*jspt - (order+1)*(order+1)*kspt;
        // Next, get evaluate Lagrange solution basis at corners
        opp_spts_to_mpts(0,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D,-1,kspt);
        opp_spts_to_mpts(1,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D,-1,kspt);
        opp_spts_to_mpts(2,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D,-1,kspt);
        opp_spts_to_mpts(3,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D,-1,kspt);

        opp_spts_to_mpts(4,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D, 1,kspt);
        opp_spts_to_mpts(5,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D, 1,kspt);
        opp_spts_to_mpts(6,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D, 1,kspt);
        opp_spts_to_mpts(7,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D, 1,kspt);

        // Last, evaluate at edge mid-nodes (to connect central flux points to edges)
        // First the bottom-face edges, then the top-face edges, then the vertical/connecting edges
        for (uint i=0; i<order+1; i++) {
          double x1 = locSpts1D[i];
          double x2 = locSpts1D[order-i];
          //                 edge                 x-location                    y-location                    z-location
          opp_spts_to_mpts(nv+i*ne+0,spt) = Lagrange(locSpts1D,x1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D,-1,kspt);
          opp_spts_to_mpts(nv+i*ne+1,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,x1,jspt) * Lagrange(locSpts1D,-1,kspt);
          opp_spts_to_mpts(nv+i*ne+2,spt) = Lagrange(locSpts1D,x2,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D,-1,kspt);
          opp_spts_to_mpts(nv+i*ne+3,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,x2,jspt) * Lagrange(locSpts1D,-1,kspt);

          opp_spts_to_mpts(nv+i*ne+4,spt) = Lagrange(locSpts1D,x1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D, 1,kspt);
          opp_spts_to_mpts(nv+i*ne+5,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,x1,jspt) * Lagrange(locSpts1D, 1,kspt);
          opp_spts_to_mpts(nv+i*ne+6,spt) = Lagrange(locSpts1D,x2,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D, 1,kspt);
          opp_spts_to_mpts(nv+i*ne+7,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,x2,jspt) * Lagrange(locSpts1D, 1,kspt);

          opp_spts_to_mpts(nv+i*ne+8,spt) = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D,x1,kspt);
          opp_spts_to_mpts(nv+i*ne+9,spt) = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,-1,jspt) * Lagrange(locSpts1D,x1,kspt);
          opp_spts_to_mpts(nv+i*ne+10,spt)= Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D,x1,kspt);
          opp_spts_to_mpts(nv+i*ne+11,spt)= Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D, 1,jspt) * Lagrange(locSpts1D,x1,kspt);
        }
      }
      break;
    }
    default:
      FatalError("Element type not yet supported.");
  }
}

matrix<double> oper::setupInterpolateSptsIpts(matrix<double> &loc_ipts)
{
  uint nIpts = loc_ipts.dim0;
  matrix<double> opp_interp(nIpts,nSpts);

  switch(eType) {
    case TRI: {
      for (uint ipt=0; ipt<nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);

        // Use the orthogonal 2D Dubiner basis for triangular elements
        for (uint spt=0; spt<nSpts; spt++)
          opp_interp(ipt,spt) = eval_dubiner_basis_2d(pt,spt,order);
      }
      break;
    }
    case QUAD: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint ipt=0; ipt<nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);
        for (uint spt=0; spt<nSpts; spt++) {
          // Structured I,J indices of current solution point
          uint ispt = spt%(nSpts/(order+1));
          uint jspt = floor(spt/(order+1));
          // 3D Tensor-Product Lagrange Interpolation
          opp_interp(ipt,spt) =  Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt);
        }
      }
      break;
    }
    case HEX: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint ipt = 0; ipt < nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);
        for (uint spt=0; spt<nSpts; spt++) {
          // Structured I,J,K indices of current solution point
          uint kspt = spt/((order+1)*(order+1));
          uint jspt = (spt-(order+1)*(order+1)*kspt)/(order+1);
          uint ispt = spt - (order+1)*jspt - (order+1)*(order+1)*kspt;
          // 3D Tensor-Product Lagrange Interpolation
          opp_interp(ipt,spt) =  Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt) * Lagrange(locSpts1D,pt.z,kspt);
        }
      }
      break;
    }
    default:
      FatalError("Element type not yet supported.");
  }

  return opp_interp;
}

void oper::getInterpWeights(double* loc_ipt, double* weights)
{
  // loc_ipt contains [x,y,z] in reference coordinates
  point pt = point(loc_ipt);

  // Note: 'weights' should be pre-size to nSpts

  switch(eType) {
    case TRI: {
      // Use the orthogonal 2D Dubiner basis for triangular elements
      for (uint spt=0; spt<nSpts; spt++)
        weights[spt] = eval_dubiner_basis_2d(pt,spt,order);
      break;
    }
    case QUAD: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint spt=0; spt<nSpts; spt++) {
        // Structured I,J indices of current solution point
        uint ispt = spt%(nSpts/(order+1));
        uint jspt = floor(spt/(order+1));
        // 3D Tensor-Product Lagrange Interpolation
        weights[spt] =  Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt);
      }
      break;
    }
    case HEX: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint spt=0; spt<nSpts; spt++) {
        // Structured I,J,K indices of current solution point
        uint kspt = spt/((order+1)*(order+1));
        uint jspt = (spt-(order+1)*(order+1)*kspt)/(order+1);
        uint ispt = spt - (order+1)*jspt - (order+1)*(order+1)*kspt;
        // 3D Tensor-Product Lagrange Interpolation
        weights[spt] =  Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt) * Lagrange(locSpts1D,pt.z,kspt);
      }
      break;
    }
    default:
      FatalError("Element type not yet supported.");
  }
}

void oper::interpolateSptsToPoints(matrix<double> &Q_spts,matrix<double> &Q_ipts, matrix<double> &loc_ipts)
{
  uint nIpts = loc_ipts.dim0;
  uint nFields = Q_spts.dim1;

  Q_ipts.setup(nIpts,nFields);
  Q_ipts.initializeToZero();

  switch(eType) {
    case TRI: {
      for (uint ipt=0; ipt<nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);

        // Use the orthogonal 2D Dubiner basis for triangular elements
        for (uint spt=0; spt<nSpts; spt++)
          for (uint field=0; field<Q_spts.dim1; field++)
            Q_ipts(ipt,field) += eval_dubiner_basis_2d(pt,spt,order);
      }
      break;
    }
    case QUAD: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint ipt=0; ipt<nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);
        for (uint spt=0; spt<nSpts; spt++) {
          // Structured I,J indices of current solution point
          uint ispt = spt%(nSpts/(order+1));
          uint jspt = floor(spt/(order+1));

          // 3D Tensor-Product Lagrange Interpolation
          for (uint field=0; field<nFields; field++)
            Q_ipts(ipt,field) += Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt);
        }
      }
      break;
    }

    case HEX: {
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint ipt = 0; ipt < nIpts; ipt++) {
        // Location of the current interpolation point
        point pt = point(loc_ipts[ipt]);
        for (uint spt=0; spt<nSpts; spt++) {
          // Structured I,J,K indices of current solution point
          uint kspt = spt/((order+1)*(order+1));
          uint jspt = (spt-(order+1)*(order+1)*kspt)/(order+1);
          uint ispt = spt - (order+1)*jspt - (order+1)*(order+1)*kspt;

          // 3D Tensor-Product Lagrange Interpolation
          for (uint field=0; field<nFields; field++)
            Q_ipts(ipt,field) += Lagrange(locSpts1D,pt.x,ispt) * Lagrange(locSpts1D,pt.y,jspt) * Lagrange(locSpts1D,pt.z,kspt);
        }
      }
      break;
    }
    default:
      FatalError("Element type not yet supported.");
  }
}

void oper::setupGradSpts(vector<point> &loc_spts)
{
  uint nSpts, spt1, spt2, dim;
  uint ispt1, jspt1, kspt1, ispt2, jspt2, kspt2;
  nSpts = loc_spts.size();

  opp_grad_spts.resize(nDims);
  for (auto& dim:opp_grad_spts) dim.setup(nSpts,nSpts);

  if (eType == TRI) {
    for (spt1=0; spt1<nSpts; spt1++) {
      for (spt2=0; spt2<nSpts; spt2++) {
        opp_grad_spts[0](spt1,spt2) = eval_dr_dubiner_basis_2d(loc_spts[spt1],spt2,order); // double-check the spt vs. spt2 in this
        opp_grad_spts[1](spt1,spt2) = eval_ds_dubiner_basis_2d(loc_spts[spt1],spt2,order);
      }
    }
  }
  else if (eType == QUAD) {
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    for (dim=0; dim<nDims; dim++) {
      for (spt1=0; spt1<nSpts; spt1++) {
        ispt1 = spt1%(nSpts/(order+1));      // col index - also = to (spt1 - (order+1)*row)
        jspt1 = floor(spt1/(order+1));       // row index
        for (spt2=0; spt2<nSpts; spt2++) {
          ispt2 = spt2%(nSpts/(order+1));
          jspt2 = floor(spt2/(order+1));
          if (dim==0) {
            opp_grad_spts[dim](spt1,spt2) = dLagrange(loc_spts_1D,loc_spts_1D[ispt1],ispt2) * Lagrange(loc_spts_1D,loc_spts_1D[jspt1],jspt2);
          } else {
            opp_grad_spts[dim](spt1,spt2) = dLagrange(loc_spts_1D,loc_spts_1D[jspt1],jspt2) * Lagrange(loc_spts_1D,loc_spts_1D[ispt1],ispt2);
          }
        }
      }
    }
  }
  else if (eType == HEX) {
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    for (dim=0; dim<nDims; dim++) {
      for (spt1=0; spt1<nSpts; spt1++) {
        kspt1 = spt1/((order+1)*(order+1));                         // col index - also = to (spt1 - (order+1)*row)
        jspt1 = (spt1-(order+1)*(order+1)*kspt1)/(order+1);         // row index
        ispt1 = spt1 - (order+1)*jspt1 - (order+1)*(order+1)*kspt1; // page index
        for (spt2=0; spt2<nSpts; spt2++) {
          kspt2 = spt2/((order+1)*(order+1));                         // col index - also = to (spt1 - (order+1)*row)
          jspt2 = (spt2-(order+1)*(order+1)*kspt2)/(order+1);         // row index
          ispt2 = spt2 - (order+1)*jspt2 - (order+1)*(order+1)*kspt2; // page index
          if (dim == 0) {
            opp_grad_spts[dim](spt1,spt2) = dLagrange(loc_spts_1D,loc_spts_1D[ispt1],ispt2) * Lagrange(loc_spts_1D,loc_spts_1D[jspt1],jspt2) * Lagrange(loc_spts_1D,loc_spts_1D[kspt1],kspt2);
          }
          else if (dim == 1) {
            opp_grad_spts[dim](spt1,spt2) = dLagrange(loc_spts_1D,loc_spts_1D[jspt1],jspt2) * Lagrange(loc_spts_1D,loc_spts_1D[ispt1],ispt2) * Lagrange(loc_spts_1D,loc_spts_1D[kspt1],kspt2);
          }
          else if (dim == 2) {
            opp_grad_spts[dim](spt1,spt2) = dLagrange(loc_spts_1D,loc_spts_1D[kspt1],kspt2) * Lagrange(loc_spts_1D,loc_spts_1D[ispt1],ispt2) * Lagrange(loc_spts_1D,loc_spts_1D[jspt1],jspt2);
          }
        }
      }
    }
  }
  else {
    FatalError("Element type not yet supported.");
  }
}


void oper::setupCorrection(vector<point> &loc_spts, vector<point> &loc_fpts)
{
  uint nSpts, nFpts;
  vector<double> loc(nDims);
  nSpts = loc_spts.size();
  nFpts = loc_fpts.size();

  opp_correction.setup(nSpts,nFpts);
  opp_correction.initializeToZero();

  if (eType == TRI) {
    // Not yet implemented
  }
  else if (eType == QUAD) {
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    for(uint spt=0; spt<nSpts; spt++) {
      for(uint fpt=0; fpt<nFpts; fpt++) {
        opp_correction(spt,fpt) = divVCJH_quad(fpt,loc_spts[spt],loc_spts_1D,params->vcjhSchemeQuad,order);
      }
    }
  }
  else if (eType == HEX) {
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    for(uint spt=0; spt<nSpts; spt++) {
      for(uint fpt=0; fpt<nFpts; fpt++) {
        opp_correction(spt,fpt) = divVCJH_hex(fpt,loc_spts[spt],loc_spts_1D,params->vcjhSchemeQuad,order);
      }
    }
  }
}

void oper::setupCorrectGradU(void)
{
  opp_correctU.resize(nDims);
  for (auto& op:opp_correctU) {
    op.setup(nSpts,nFpts);
    op.initializeToZero();
  }

  if (eType == TRI) {
    // Not yet implemented
  }
  else if (eType == QUAD) {

    vector<double> tNorm(nDims);

    for(uint fpt=0; fpt<nFpts; fpt++) {

      uint iFace = floor(fpt / (order+1));
      switch(iFace) {
      case(0):
        tNorm[0] = 0;
        tNorm[1] = -1;
        break;
      case(1):
        tNorm[0] = 1;
        tNorm[1] = 0;
        break;
      case(2):
        tNorm[0] = 0;
        tNorm[1] = 1;
        break;
      case(3):
        tNorm[0] = -1;
        tNorm[1] = 0;
        break;
      }

      for(uint spt=0; spt<nSpts; spt++) {
        for(uint dim=0; dim<nDims; dim++) {
          opp_correctU[dim](spt,fpt) = opp_correction(spt,fpt) * tNorm[dim];
        }
      }

    }
  }
  else if (eType == HEX) {
    vector<double> tNorm(nDims);

    for(uint fpt=0; fpt<nFpts; fpt++) {

      uint iFace = floor(fpt / ((order+1)*(order+1)));
      switch(iFace) {
        case 0:
          tNorm[0] =  0;
          tNorm[1] =  0;
          tNorm[2] = -1;
          break;
        case 1:
          tNorm[0] =  0;
          tNorm[1] =  0;
          tNorm[2] =  1;
          break;
        case 2:
          tNorm[0] = -1;
          tNorm[1] =  0;
          tNorm[2] =  0;
          break;
        case 3:
          tNorm[0] =  1;
          tNorm[1] =  0;
          tNorm[2] =  0;
          break;
        case 4:
          tNorm[0] =  0;
          tNorm[1] = -1;
          tNorm[2] =  0;
          break;
        case 5:
          tNorm[0] =  0;
          tNorm[1] =  1;
          tNorm[2] =  0;
          break;
      }

      for(uint spt=0; spt<nSpts; spt++) {
        for(uint dim=0; dim<nDims; dim++) {
          opp_correctU[dim](spt,fpt) = opp_correction(spt,fpt) * tNorm[dim];
        }
      }
    }
  }
}

// Setup Vandermonde Matrices
void oper::setupVandermonde(vector<point> &loc_spts)
{
  if (eType == TRI) {
    // Not yet implemented
  }
  else if (eType == QUAD) {
    vandermonde1D.setup(order+1,order+1);
    vandermonde2D.setup(nSpts,nSpts);
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    vector<double> loc(nDims);

    // Create 1D Vandermonde matrix
    for (uint i=0; i<order+1; i++)
      for (uint j=0; j<order+1; j++)
        vandermonde1D(i,j) = Legendre(loc_spts_1D[i],j);

    // Store its inverse
    inv_vandermonde1D = vandermonde1D.invertMatrix();

    // Create 2D Vandermonde matrix
    for (uint i=0; i<nSpts; i++) {
      loc[0] = loc_spts[i].x;
      loc[1] = loc_spts[i].y;

      for (uint j=0; j<nSpts; j++)
        vandermonde2D(i,j) = Legendre2D_hierarchical(j,loc,order);
    }

    // Store its inverse
    inv_vandermonde2D = vandermonde2D.invertMatrix();
  }
}

// Set the 1D concentration matrix based on 1D-loc_spts
void oper::setupSensingMatrix(void)
{
  if (eType == TRI) {
    // Not yet implemented
  }
  else if (eType == QUAD) {
    int concen_type = 1;    //Considering getting this from input file
    matrix<double> concentration_factor(order+1,1);
    matrix<double> grad_vandermonde;
    grad_vandermonde.setup(order+1,order+1);
    matrix<double> concentrationMatrix(order+1,order+1);

    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);

    // create the vandermonde matrix
    for(uint i=0; i<order+1; i++)
      for (uint j=0; j<order+1; j++)
        grad_vandermonde(i,j) = dLegendre(loc_spts_1D[i],j);

    // create concentration factor array
    for(uint j=0; j <order+1; j++){
      if(concen_type == 0){ // exponential
        if(j==0)
          concentration_factor(j) = 0;
        else
          concentration_factor(j) = exp(1/(6*j*(j+1)));
      }
      else if(concen_type == 1) // linear
        concentration_factor(j) = 1;

      else
        cout<<"Concentration factor not setup"<<endl;
    }

    // Prepare concentration matrix as in paper
    for(uint i=0; i<order+1; i++)
      for(uint j=0; j<order+1; j++)
        concentrationMatrix(i,j) = (3.1415/(order+1))*concentration_factor(j)*sqrt(1 - loc_spts_1D[i]*loc_spts_1D[i])*grad_vandermonde(i,j);

    concentrationMatrix.timesMatrix(inv_vandermonde1D,sensingMatrix);
  }
}

void oper::setupFilterMatrix(void)
{
  if (eType == TRI) {
    // Not yet implemented
  }
  else if (eType == QUAD) {
    matrix<double> filterWeights(nSpts,nSpts);
    matrix<double> tempMatrix(nSpts,nSpts);
    filterMatrix.setup(nSpts,nSpts);
    double exponent = 1;     // Governs filter strength - Take as input?

    // create the filter weights matrix as a diagonal matrix
    for (uint i=0; i<nSpts; i++)
      filterWeights(i,i) = exponential_filter(i,order,exponent);

    filterWeights.print();

    // Filter matrix is SigmaMatrix*V_inverse (so it can directly
    // operate on nodal solution)
    filterWeights.timesMatrix(inv_vandermonde2D,tempMatrix);
    vandermonde2D.timesMatrix(tempMatrix,filterMatrix);
  }
}

void oper::applyGradSpts(matrix<double> &U_spts, vector<matrix<double> > &dU_spts)
{
  for (uint dim=0; dim<nDims; dim++)
    opp_grad_spts[dim].timesMatrix(U_spts,dU_spts[dim]);
}

void oper::applyGradFSpts(vector<matrix<double>> &F_spts, vector<vector<matrix<double>>> &dF_spts)
{
  // Note: dim1 is flux direction, dim2 is derivative direction
  for (uint dim1=0; dim1<nDims; dim1++)
    for (uint dim2=0; dim2<dF_spts.size(); dim2++)
      opp_grad_spts[dim2].timesMatrix(F_spts[dim1],dF_spts[dim2][dim1]);
}


void oper::applyDivFSpts(vector<matrix<double>> &F_spts, matrix<double> &divF_spts)
{
  divF_spts.initializeToZero();
  for (uint dim=0; dim<nDims; dim++)
      opp_grad_spts[dim].timesMatrixPlus(F_spts[dim],divF_spts);
}


void oper::applySptsFpts(matrix<double> &U_spts, matrix<double> &U_fpts)
{
  opp_spts_to_fpts.timesMatrix(U_spts,U_fpts);
}

void oper::applySptsMpts(matrix<double> &U_spts, matrix<double> &U_mpts)
{
  opp_spts_to_mpts.timesMatrix(U_spts,U_mpts);
}

void oper::applyExtrapolateFn(vector<matrix<double>> &F_spts, matrix<double> &tnorm_fpts, matrix<double> &Fn_fpts)
{
  uint nFpts = tnorm_fpts.getDim0();
  matrix<double> tempFn(nFpts,nDims);
  tempFn.initializeToZero();
  Fn_fpts.initializeToZero();

  for (uint dim=0; dim<nDims; dim++) {
    opp_spts_to_fpts.timesMatrix(F_spts[dim],tempFn);
    for (uint fpt=0; fpt<nFpts; fpt++)
      for (uint i=0; i<nFields; i++)
        Fn_fpts(fpt,i) += tempFn(fpt,i)*tnorm_fpts(fpt,dim);
  }
}

void oper::applyExtrapolateFn(vector<matrix<double>> &F_spts, matrix<double> &norm_fpts, matrix<double> &Fn_fpts, vector<double>& dA_fpts)
{
  uint nFpts = norm_fpts.getDim0();
  matrix<double> tempFn(nFpts,nDims);
  tempFn.initializeToZero();
  Fn_fpts.initializeToZero();

  for (uint dim=0; dim<nDims; dim++) {
    opp_spts_to_fpts.timesMatrix(F_spts[dim],tempFn);
    for (uint fpt=0; fpt<nFpts; fpt++)
      for (uint i=0; i<nFields; i++)
        Fn_fpts(fpt,i) += tempFn(fpt,i)*norm_fpts(fpt,dim)*dA_fpts[fpt];
  }
}

void oper::applyCorrectDivF(matrix<double> &dFn_fpts, matrix<double> &divF_spts)
{
  opp_correction.timesMatrixPlus(dFn_fpts,divF_spts);
}

void oper::applyCorrectGradU(matrix<double> &dUc_fpts, vector<matrix<double>> &dU_spts)
{
  for (uint dim=0; dim<nDims; dim++)
    opp_correctU[dim].timesMatrixPlus(dUc_fpts,dU_spts[dim]);
}

void oper::calcAvgU(matrix<double> &U_spts, vector<double> &detJ_spts, vector<double> &Uavg)
{
  auto weights = Geo->getQptWeights(order);

  Uavg.assign(nFields,0);
  double vol = 0;
  for (uint spt=0; spt<nSpts; spt++) {
    for (uint i=0; i<nFields; i++) {
      Uavg[i] += U_spts(spt,i)*weights[spt]*detJ_spts[spt];
    }
    vol += weights[spt]*detJ_spts[spt];
  }

  for (auto &i:Uavg) i/= vol;
}


const matrix<double> &oper::get_oper_div_spts()
{
  return opp_div_spts;
}


double oper::divVCJH_quad(int in_fpt, point& loc, vector<double>& loc_1d_spts, uint vcjh, uint order)
{
  uint i,j;
  double eta;
  double div_vcjh_basis = 0;

  if (vcjh == DG)
    eta = 0.; // HiFiLES: run_input.eta_quad;
  else
    eta = compute_eta(vcjh,order);

  i = in_fpt / (order+1);      // Face upon which the flux point lies [0,1,2, or 3]
  j = in_fpt - (order+1)*i;    // Face-local index of flux point [0 to n_fpts_per_face-1]

  if(i==0)       // Bottom
    div_vcjh_basis = -Lagrange(loc_1d_spts,loc[0],j) * dVCJH_1d(loc[1],0,order,eta); // was -'ve
  else if(i==1)  // Right
    div_vcjh_basis =  Lagrange(loc_1d_spts,loc[1],j) * dVCJH_1d(loc[0],1,order,eta);
  else if(i==2)  // Top
    div_vcjh_basis =  Lagrange(loc_1d_spts,loc[0],order-j) * dVCJH_1d(loc[1],1,order,eta);
  else if(i==3)  // Left
    div_vcjh_basis = -Lagrange(loc_1d_spts,loc[1],order-j) * dVCJH_1d(loc[0],0,order,eta); // was -'ve


  return div_vcjh_basis;
}

double oper::divVCJH_hex(int in_fpt, point& loc, vector<double>& loc_1d_spts, uint vcjh, uint order)
{
  double eta;
  double div_vcjh_basis = 0;

  if (vcjh == DG)
    eta = 0.; // HiFiLES: run_input.eta_quad;
  else
    eta = compute_eta(vcjh,order);

  uint P1  = order+1;
  uint P12 = P1*P1;
  uint f = floor(in_fpt / P12);     // Face upon which the flux point lies [0 to 5]
  uint j = (in_fpt - P12*f) / P1;   // Face-local index of flux point [0 to n_fpts_per_face-1]
  uint i =  in_fpt - P12*f  - P1*j; // Face-local index of flux point [0 to n_fpts_per_face-1]

  if(f==0)       // z-plane, 'L'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.x,i)       * Lagrange(loc_1d_spts,loc.y,j) * -dVCJH_1d(loc.z,0,order,eta);
  else if(f==1)  // z-plane, 'R'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.x,order-i) * Lagrange(loc_1d_spts,loc.y,j) *  dVCJH_1d(loc.z,1,order,eta);
  else if(f==2)  // x-plane, 'L'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.y,i)       * Lagrange(loc_1d_spts,loc.z,j) * -dVCJH_1d(loc.x,0,order,eta);
  else if(f==3)  // x-plane, 'R'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.y,order-i) * Lagrange(loc_1d_spts,loc.z,j) *  dVCJH_1d(loc.x,1,order,eta);
  else if(f==4)  // y-plane, 'L'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.x,order-i) * Lagrange(loc_1d_spts,loc.z,j) * -dVCJH_1d(loc.y,0,order,eta);
  else if(f==5)  // y-plane, 'R'
    div_vcjh_basis = Lagrange(loc_1d_spts,loc.x,i)       * Lagrange(loc_1d_spts,loc.z,j) *  dVCJH_1d(loc.y,1,order,eta);

  return div_vcjh_basis;
}

// Method to capture shock in the element
double oper::shockCaptureInEle(matrix<double> &U_spts, double threshold)
{
  double sensor = 0;
  if(eType == TRI)
    cout << "Shock capturing not implemented yet for Dubiner-basis triangles" << endl;
  else if (eType == QUAD) {

    // Sensing Part
    int p = 3;  // Exponent of concentration method
    matrix<double> Rho(order+1,order+1);
    vector<double> rho = U_spts.getCol(0);
    vector<double> rho_x, rho_y, uEx, uEy;
    std::vector<double>::iterator maxuEx, maxuEy;
    double newmax, currmax = 0;

    // Put the density into a 2D matrix- to access as row,col
    Rho.vecToMatrixResize(rho);

    // Sense along the X and Y-slices
    for(uint i=0; i<order+1; i++) {
      rho_x = Rho.getRow(i);
      rho_y = Rho.getCol(i);
      sensingMatrix.timesVector(rho_x,uEx);
      sensingMatrix.timesVector(rho_x,uEy);
      maxuEx = max_element(uEx.begin(), uEx.end(), abs_compare);
      maxuEy = max_element(uEx.begin(), uEx.end(), abs_compare);
      newmax = max(abs(*maxuEx),abs(*maxuEy));
      currmax = max(currmax,newmax);
    }

    sensor = pow(currmax,p)*pow(order+1,p/2);

    //Filtering Part
    if(sensor >= threshold) {
      //cout<<"It is filtering!"<<endl;
      matrix<double> tempMatrix(U_spts.getDim0(),U_spts.getDim1());
      filterMatrix.timesMatrix(U_spts,tempMatrix);
      U_spts = tempMatrix;
    }
  }
  return sensor;
}
