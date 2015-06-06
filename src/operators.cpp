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

  // Operators needed for Shock capturing
  if (params->scFlag) {
    setupVandermonde(loc_spts);

    setupSensingMatrix();

    setupFilterMatrix();
  }
}

void oper::setupExtrapolateSptsFpts(vector<point> &loc_fpts)
{
  uint spt, fpt, ispt, jspt;
  opp_spts_to_fpts.setup(nFpts,nSpts);

  for (fpt=0; fpt<nFpts; fpt++) {
    for (spt=0; spt<nSpts; spt++) {
      switch(eType) {
        case(TRI):
          opp_spts_to_fpts[fpt][spt] = eval_dubiner_basis_2d(loc_fpts[fpt],spt,order);
          break;
        case(QUAD): {
          vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
          // First, get the i an j ID of the spt
          ispt = spt%(nSpts/(order+1));
          jspt = floor(spt/(order+1));
          opp_spts_to_fpts[fpt][spt] = Lagrange(locSpts1D,loc_fpts[fpt].x,ispt) * Lagrange(locSpts1D,loc_fpts[fpt].y,jspt);
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
    case(TRI): {
      opp_spts_to_mpts.setup(3,nSpts);
      point vert1, vert2, vert3;
      vert1.x = -1;  vert1.y = -1;
      vert2.x =  1;  vert2.y = -1;
      vert3.x = -1;  vert3.y =  1;
      for (uint spt=0; spt<nSpts; spt++) {
        opp_spts_to_mpts[0][spt] = eval_dubiner_basis_2d(vert1,spt,order);
        opp_spts_to_mpts[1][spt] = eval_dubiner_basis_2d(vert2,spt,order);
        opp_spts_to_mpts[2][spt] = eval_dubiner_basis_2d(vert3,spt,order);
      }
      break;
    }
    case(QUAD): {
      uint ispt, jspt;
      opp_spts_to_mpts.setup(4,nSpts);
      vector<double> locSpts1D = Geo->getPts1D(params->sptsTypeQuad,order);
      for (uint spt=0; spt<nSpts; spt++) {
        // First, get the i an j ID of the spt
        ispt = spt%(nSpts/(order+1));
        jspt = floor(spt/(order+1));
        // Next, get evaluate Lagrange solution basis at corners
        opp_spts_to_mpts[0][spt] = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D,-1,jspt);
        opp_spts_to_mpts[1][spt] = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D,-1,jspt);
        opp_spts_to_mpts[2][spt] = Lagrange(locSpts1D, 1,ispt) * Lagrange(locSpts1D, 1,jspt);
        opp_spts_to_mpts[3][spt] = Lagrange(locSpts1D,-1,ispt) * Lagrange(locSpts1D, 1,jspt);
      }
      break;
    }
    default:
      FatalError("Element type not yet supported.");
  }
}

void oper::setupInterpolate(vector<point> &pts_from, vector<point> &pts_to, matrix<double> &opp_interp)
{
  uint ptA, ptB, nPtsFrom, nPtsTo, pt, iptA, jptA;
  nPtsFrom = pts_from.size();
  nPtsTo = pts_to.size();

  opp_interp.setup(nPtsFrom,nPtsTo);

  // Get 1D locations of points for arbitrary, potentially anisotropic tensor-product elements
  vector<double> locPts1Dx, locPts1Dy;
  if (eType==QUAD) {
    locPts1Dx.resize(order+1);
    locPts1Dy.resize(order+1);
    for (pt=0; pt<(uint)order+1; pt++) {
      locPts1Dx[pt] = pts_from[pt].x;
      locPts1Dy[pt] = pts_from[pt*(order+1)].y;
    }
  }

  for (ptB=0; ptB<nPtsTo; ptB++) {
    for (ptA=0; ptA<nPtsFrom; ptA++) {
      switch(eType) {
        case(TRI):
          opp_interp[ptB][ptA] = eval_dubiner_basis_2d(pts_to[ptB],ptA,order);
          break;
        case(QUAD):
          // First, get the i and j ID of the pt [tensor-product element]
          iptA = ptA%(nPtsFrom/(order+1));
          jptA = floor(ptA/(order+1));
          opp_interp[ptB][ptA] = Lagrange(locPts1Dx,pts_to[ptB].x,iptA) * Lagrange(locPts1Dy,pts_to[ptB].y,jptA);
          break;
        default:
          FatalError("Element type not yet supported.");
      }
    }
  }
}

void oper::setupGradSpts(vector<point> &loc_spts)
{
  uint nSpts, spt1, spt2, dim;
  uint xid1, yid1, xid2, yid2;
  nSpts = loc_spts.size();

  opp_grad_spts.resize(nDims);
  for (auto& dim:opp_grad_spts) dim.setup(nSpts,nSpts);

  if (eType == TRI) {
    for (spt1=0; spt1<nSpts; spt1++) {
      for (spt2=0; spt2<nSpts; spt2++) {
        opp_grad_spts[0][spt1][spt2] = eval_dr_dubiner_basis_2d(loc_spts[spt1],spt2,order); // double-check the spt vs. spt2 in this
        opp_grad_spts[1][spt1][spt2] = eval_ds_dubiner_basis_2d(loc_spts[spt1],spt2,order);
      }
    }
  }
  else if (eType == QUAD) {
    vector<double> loc_spts_1D = Geo->getPts1D(params->sptsTypeQuad,order);
    for (dim=0; dim<nDims; dim++) {
      for (spt1=0; spt1<nSpts; spt1++) {
        xid1 = spt1%(nSpts/(order+1));      // col index - also = to (spt1 - (order+1)*row)
        yid1 = floor(spt1/(order+1));       // row index
        for (spt2=0; spt2<nSpts; spt2++) {
          xid2 = spt2%(nSpts/(order+1));
          yid2 = floor(spt2/(order+1));
          if (dim==0) {
            opp_grad_spts[dim][spt1][spt2] = dLagrange(loc_spts_1D,loc_spts_1D[xid1],xid2) * Lagrange(loc_spts_1D,loc_spts_1D[yid1],yid2);
          } else {
            opp_grad_spts[dim][spt1][spt2] = dLagrange(loc_spts_1D,loc_spts_1D[yid1],yid2) * Lagrange(loc_spts_1D,loc_spts_1D[xid1],xid2);
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
  uint nSpts, nFpts, spt, fpt, dim;
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
    for(spt=0; spt<nSpts; spt++) {
      for(dim=0; dim<nDims; dim++){
        loc[dim] = loc_spts[spt][dim];
      }

      for(fpt=0; fpt<nFpts; fpt++) {
        opp_correction[spt][fpt] = divVCJH_quad(fpt,loc,loc_spts_1D,params->vcjhSchemeQuad,order);
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
        Fn_fpts[fpt][i] += tempFn[fpt][i]*tnorm_fpts[fpt][dim];
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
        Fn_fpts[fpt][i] += tempFn[fpt][i]*norm_fpts[fpt][dim]*dA_fpts[fpt];
  }
}

void oper::applyCorrectDivF(matrix<double> &dFn_fpts, matrix<double> &divF_spts)
{
  opp_correction.timesMatrixPlus(dFn_fpts,divF_spts);
}


const matrix<double> &oper::get_oper_div_spts()
{
  return opp_div_spts;
}


double oper::divVCJH_quad(int in_fpt, vector<double>& loc, vector<double>& loc_1d_spts, uint vcjh, uint order)
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

  if(i==0)
    div_vcjh_basis = -Lagrange(loc_1d_spts,loc[0],j) * dVCJH_1d(loc[1],0,order,eta); // was -'ve
  else if(i==1)
    div_vcjh_basis =  Lagrange(loc_1d_spts,loc[1],j) * dVCJH_1d(loc[0],1,order,eta);
  else if(i==2)
    div_vcjh_basis =  Lagrange(loc_1d_spts,loc[0],order-j) * dVCJH_1d(loc[1],1,order,eta);
  else if(i==3)
    div_vcjh_basis =-Lagrange(loc_1d_spts,loc[1],order-j) * dVCJH_1d(loc[0],0,order,eta); // was -'ve


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
