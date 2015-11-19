/*!
 * \file funcs.cpp
 * \brief Miscellaneous helper functions
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
#include "funcs.hpp"

vector<double> solveCholesky(matrix<double> A, vector<double> b)
{
  int m = A.getDim0();
  int n = A.getDim1();
  if (m!=n) FatalError("Cannot use Cholesky on non-square matrix.");

  // Get the Cholesky factorization of A [A = G*G^T]
  for (int j=0; j<n; j++) {
    if (j>0) {
      matrix<double> a = A.slice({{j,n-1}},{{0,j-1}});
      vector<double> a1 = a.getRow(0);
      a1 = a*a1;
      for (int i=0; i<n-j; i++) A(j+i,j) -= a1[i];
    }

    if (A(j,j)<0) {
      A.print();
      FatalError("Negative factor in Cholesky!");
    }

    double ajj = sqrt(A(j,j));
    for (int i=j; i<n; i++)
      A(i,j) /= ajj;
  }

  // Lower-Triangular Solve [G*y = b]
  vector<double> y(n);
  for (int i=0; i<n; i++) {
    y[i] = b[i];
    for (int j=0; j<i; j++) {
      y[i] -= A(i,j)*y[j];
    }
    y[i] /= A(i,i);
  }

  // Upper-Triangular Solve [G^T*x = y]
  vector<double> x(n);
  for (int i=n-1; i>=0; i--) {
    x[i] = y[i];
    for (int j=i+1; j<n; j++) {
      x[i] -= A(j,i)*x[j];
    }
    x[i] /= A(i,i);
  }

  return x;
}

matrix<double> solveCholesky(matrix<double> A, matrix<double> &B)
{
  double eps = 1e-12;
  int m = A.getDim0();
  int n = A.getDim1();
  int p = B.getDim1();
  if (m!=n) FatalError("Cannot use Cholesky on non-square matrix.");

  // Get the Cholesky factorization of A [A = G*G^T]
  for (int j=0; j<n; j++) {
    if (j>0) {
      matrix<double> a = A.slice({{j,n-1}},{{0,j-1}});
      vector<double> a1 = a.getRow(0);
      a1 = a*a1;
      for (int i=0; i<n-j; i++) A(j+i,j) -= a1[i];
    }

    if (A(j,j)<0) {
      if (std::abs(A(j,j) < eps)) {
        A(j,j) = eps;
      } else {
        _(A(j,j));
        FatalError("Negative factor in Cholesky!");
      }
    }

    double ajj = sqrt(A(j,j));
    for (int i=j; i<n; i++) {
      A(i,j) /= ajj;
    }
  }

  // Lower-Triangular Solve [G*Y = B]
  matrix<double> y(n,p);
  for (int k=0; k<p; k++) {
    for (int i=0; i<n; i++) {
      y(i,k) = B(i,k);
      for (int j=0; j<i; j++) {
        y(i,k) -= A(i,j)*y(j,k);
      }
      y(i,k) /= A(i,i);
    }
  }

  // Upper-Triangular Solve [G^T*X = Y]
  matrix<double> x(n,p);
  for (int k=0; k<p; k++) {
    for (int i=n-1; i>=0; i--) {
      x(i,k) = y(i,k);
      for (int j=i+1; j<n; j++) {
        x(i,k) -= A(j,i)*x(j,k);
      }
      x(i,k) /= A(i,i);
    }
  }

  return x;
}

void shape_quad(const point &in_rs, vector<double> &out_shape, int nNodes)
{
  out_shape.resize(nNodes);
  shape_quad(in_rs, out_shape.data(), nNodes);
}

void shape_quad(const point &in_rs, double* out_shape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;
  switch(nNodes) {
  case 4:
    out_shape[0] = 0.25*(1-xi)*(1-eta);
    out_shape[1] = 0.25*(1+xi)*(1-eta);
    out_shape[2] = 0.25*(1+xi)*(1+eta);
    out_shape[3] = 0.25*(1-xi)*(1+eta);
    break;

  case 8:
    out_shape[0] = -0.25*(1-xi)*(1-eta)*(1+eta+xi);
    out_shape[1] = -0.25*(1+xi)*(1-eta)*(1+eta-xi);
    out_shape[2] = -0.25*(1+xi)*(1+eta)*(1-eta-xi);
    out_shape[3] = -0.25*(1-xi)*(1+eta)*(1-eta+xi);
    out_shape[4] = 0.5*(1-xi)*(1+ xi)*(1-eta);
    out_shape[5] = 0.5*(1+xi)*(1+eta)*(1-eta);
    out_shape[6] = 0.5*(1-xi)*(1+ xi)*(1+eta);
    out_shape[7] = 0.5*(1-xi)*(1+eta)*(1-eta);
    break;
  }
}

void shape_hex(const point &in_rst, vector<double> &out_shape, int nNodes)
{
  out_shape.resize(nNodes);
  shape_hex(in_rst, out_shape.data(), nNodes);
}

void shape_hex(const point &in_rst, double* out_shape, int nNodes)
{
  double xi  = in_rst.x;
  double eta = in_rst.y;
  double mu = in_rst.z;
  switch(nNodes) {
    case 8:
      out_shape[0] = 0.125*(1-xi)*(1-eta)*(1-mu);
      out_shape[1] = 0.125*(1+xi)*(1-eta)*(1-mu);
      out_shape[2] = 0.125*(1+xi)*(1+eta)*(1-mu);
      out_shape[3] = 0.125*(1-xi)*(1+eta)*(1-mu);

      out_shape[4] = 0.125*(1-xi)*(1-eta)*(1+mu);
      out_shape[5] = 0.125*(1+xi)*(1-eta)*(1+mu);
      out_shape[6] = 0.125*(1+xi)*(1+eta)*(1+mu);
      out_shape[7] = 0.125*(1-xi)*(1+eta)*(1+mu);
      break;

    case 20: {
      double XI[8]  = {-1,1,1,-1,-1,1,1,-1};
      double ETA[8] = {-1,-1,1,1,-1,-1,1,1};
      double MU[8]  = {-1,-1,-1,-1,1,1,1,1};
      // Corner nodes
      for (int i=0; i<8; i++) {
        out_shape[i] = .125*(1+xi*XI[i])*(1+eta*ETA[i])*(1+mu*MU[i])*(xi*XI[i]+eta*ETA[i]+mu*MU[i]-2);
      }
      // Edge nodes, xi = 0
      out_shape[8]  = .25*(1-xi*xi)*(1-eta)*(1-mu);
      out_shape[10] = .25*(1-xi*xi)*(1+eta)*(1-mu);
      out_shape[16] = .25*(1-xi*xi)*(1-eta)*(1+mu);
      out_shape[18] = .25*(1-xi*xi)*(1+eta)*(1+mu);
      // Edge nodes, eta = 0
      out_shape[9]  = .25*(1-eta*eta)*(1+xi)*(1-mu);
      out_shape[11] = .25*(1-eta*eta)*(1-xi)*(1-mu);
      out_shape[17] = .25*(1-eta*eta)*(1+xi)*(1+mu);
      out_shape[19] = .25*(1-eta*eta)*(1-xi)*(1+mu);
      // Edge Nodes, mu = 0
      out_shape[12] = .25*(1-mu*mu)*(1-xi)*(1-eta);
      out_shape[13] = .25*(1-mu*mu)*(1+xi)*(1-eta);
      out_shape[14] = .25*(1-mu*mu)*(1+xi)*(1+eta);
      out_shape[15] = .25*(1-mu*mu)*(1-xi)*(1+eta);
      break;
    }
  }
}

void dshape_quad(const point &in_rs, matrix<double> &out_dshape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;
  out_dshape.setup(nNodes,2);

  switch(nNodes) {
  case 4:
    out_dshape(0,0) = -0.25*(1-eta);
    out_dshape(1,0) =  0.25*(1-eta);
    out_dshape(2,0) =  0.25*(1+eta);
    out_dshape(3,0) = -0.25*(1+eta);

    out_dshape(0,1) = -0.25*(1-xi);
    out_dshape(1,1) = -0.25*(1+xi);
    out_dshape(2,1) =  0.25*(1+xi);
    out_dshape(3,1) =  0.25*(1-xi);
    break;

  case 8:
    out_dshape(0,0) = -0.25*(-1+eta)*(2*xi+eta);
    out_dshape(1,0) =  0.25*(-1+eta)*(eta - 2*xi);
    out_dshape(2,0) =  0.25*( 1+eta)*(2*xi+eta);
    out_dshape(3,0) = -0.25*( 1+eta)*(eta-2*xi);
    out_dshape(4,0) =    xi*(-1+eta);
    out_dshape(5,0) = -0.5 *( 1+eta)*(-1+eta);
    out_dshape(6,0) =   -xi*( 1+eta);
    out_dshape(7,0) =  0.5 *( 1+eta)*(-1+eta);

    out_dshape(0,1) = -0.25*(-1+xi)*(2*eta+xi);
    out_dshape(1,1) =  0.25*( 1+xi)*(2*eta - xi);
    out_dshape(2,1) =  0.25*( 1+xi)*(2*eta+xi);
    out_dshape(3,1) = -0.25*(-1+xi)*(2*eta-xi);
    out_dshape(4,1) =  0.5 *( 1+xi)*(-1+xi);
    out_dshape(5,1) =  -eta*( 1+xi);
    out_dshape(6,1) = -0.5 *( 1+xi)*(-1+xi);
    out_dshape(7,1) =   eta*(-1+xi);
    break;
  }
}

void dshape_hex(const point &in_rst, matrix<double> &out_dshape, int nNodes)
{
  double xi  = in_rst.x;
  double eta = in_rst.y;
  double mu = in_rst.z;
  out_dshape.setup(nNodes,3);

  switch(nNodes) {
    case 8:
      out_dshape(0,0) = -0.125*(1-eta)*(1-mu);
      out_dshape(1,0) =  0.125*(1-eta)*(1-mu);
      out_dshape(2,0) =  0.125*(1+eta)*(1-mu);
      out_dshape(3,0) = -0.125*(1+eta)*(1-mu);

      out_dshape(4,0) = -0.125*(1-eta)*(1+mu);
      out_dshape(5,0) =  0.125*(1-eta)*(1+mu);
      out_dshape(6,0) =  0.125*(1+eta)*(1+mu);
      out_dshape(7,0) = -0.125*(1+eta)*(1+mu);

      out_dshape(0,1) = -0.125*(1-xi)*(1-mu);
      out_dshape(1,1) = -0.125*(1+xi)*(1-mu);
      out_dshape(2,1) =  0.125*(1+xi)*(1-mu);
      out_dshape(3,1) =  0.125*(1-xi)*(1-mu);

      out_dshape(4,1) = -0.125*(1-xi)*(1+mu);
      out_dshape(5,1) = -0.125*(1+xi)*(1+mu);
      out_dshape(6,1) =  0.125*(1+xi)*(1+mu);
      out_dshape(7,1) =  0.125*(1-xi)*(1+mu);

      out_dshape(0,2) = -0.125*(1-xi)*(1-eta);
      out_dshape(1,2) = -0.125*(1+xi)*(1-eta);
      out_dshape(2,2) = -0.125*(1+xi)*(1+eta);
      out_dshape(3,2) = -0.125*(1-xi)*(1+eta);

      out_dshape(4,2) =  0.125*(1-xi)*(1-eta);
      out_dshape(5,2) =  0.125*(1+xi)*(1-eta);
      out_dshape(6,2) =  0.125*(1+xi)*(1+eta);
      out_dshape(7,2) =  0.125*(1-xi)*(1+eta);
      break;
    case 20: {
      double XI[8]  = {-1,1,1,-1,-1,1,1,-1};
      double ETA[8] = {-1,-1,1,1,-1,-1,1,1};
      double MU[8]  = {-1,-1,-1,-1,1,1,1,1};
      // Corner Nodes
      for (int i=0; i<8; i++) {
        out_dshape(i,0) = .125*XI[i] *(1+eta*ETA[i])*(1 + mu*MU[i])*(2*xi*XI[i] +   eta*ETA[i] +   mu*MU[i]-1);
        out_dshape(i,1) = .125*ETA[i]*(1 + xi*XI[i])*(1 + mu*MU[i])*(  xi*XI[i] + 2*eta*ETA[i] +   mu*MU[i]-1);
        out_dshape(i,2) = .125*MU[i] *(1 + xi*XI[i])*(1+eta*ETA[i])*(  xi*XI[i] +   eta*ETA[i] + 2*mu*MU[i]-1);
      }
      // Edge Nodes, xi = 0
      out_dshape( 8,0) = -.5*xi*(1-eta)*(1-mu);  out_dshape( 8,1) = -.25*(1-xi*xi)*(1-mu);  out_dshape( 8,2) = -.25*(1-xi*xi)*(1-eta);
      out_dshape(10,0) = -.5*xi*(1+eta)*(1-mu);  out_dshape(10,1) =  .25*(1-xi*xi)*(1-mu);  out_dshape(10,2) = -.25*(1-xi*xi)*(1+eta);
      out_dshape(16,0) = -.5*xi*(1-eta)*(1+mu);  out_dshape(16,1) = -.25*(1-xi*xi)*(1+mu);  out_dshape(16,2) =  .25*(1-xi*xi)*(1-eta);
      out_dshape(18,0) = -.5*xi*(1+eta)*(1+mu);  out_dshape(18,1) =  .25*(1-xi*xi)*(1+mu);  out_dshape(18,2) =  .25*(1-xi*xi)*(1+eta);
      // Edge Nodes, eta = 0
      out_dshape( 9,1) = -.5*eta*(1+xi)*(1-mu);  out_dshape( 9,0) =  .25*(1-eta*eta)*(1-mu);  out_dshape( 9,2) = -.25*(1-eta*eta)*(1+xi);
      out_dshape(11,1) = -.5*eta*(1-xi)*(1-mu);  out_dshape(11,0) = -.25*(1-eta*eta)*(1-mu);  out_dshape(11,2) = -.25*(1-eta*eta)*(1-xi);
      out_dshape(17,1) = -.5*eta*(1+xi)*(1+mu);  out_dshape(17,0) =  .25*(1-eta*eta)*(1+mu);  out_dshape(17,2) =  .25*(1-eta*eta)*(1+xi);
      out_dshape(19,1) = -.5*eta*(1-xi)*(1+mu);  out_dshape(19,0) = -.25*(1-eta*eta)*(1+mu);  out_dshape(19,2) =  .25*(1-eta*eta)*(1-xi);
      // Edge Nodes, mu = 0;
      out_dshape(12,2) = -.5*mu*(1-xi)*(1-eta);  out_dshape(12,0) = -.25*(1-mu*mu)*(1-eta);  out_dshape(12,1) = -.25*(1-mu*mu)*(1-xi);
      out_dshape(13,2) = -.5*mu*(1+xi)*(1-eta);  out_dshape(13,0) =  .25*(1-mu*mu)*(1-eta);  out_dshape(13,1) = -.25*(1-mu*mu)*(1+xi);
      out_dshape(14,2) = -.5*mu*(1+xi)*(1+eta);  out_dshape(14,0) =  .25*(1-mu*mu)*(1+eta);  out_dshape(14,1) =  .25*(1-mu*mu)*(1+xi);
      out_dshape(15,2) = -.5*mu*(1-xi)*(1+eta);  out_dshape(15,0) = -.25*(1-mu*mu)*(1+eta);  out_dshape(15,1) =  .25*(1-mu*mu)*(1-xi);
      break;
    }
  }
}

void ddshape_quad(const point &in_rs, Array<double,3> &out_dshape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;
  out_dshape.setup(nNodes,2,2);

  switch(nNodes) {
  case 4:
    out_dshape(0,0,0) = 0;
    out_dshape(1,0,0) = 0;
    out_dshape(2,0,0) = 0;
    out_dshape(3,0,0) = 0;

    out_dshape(0,0,1) =  0.25;
    out_dshape(1,0,1) = -0.25;
    out_dshape(2,0,1) =  0.25;
    out_dshape(3,0,1) = -0.25;

    out_dshape(0,1,0) =  0.25;
    out_dshape(1,1,0) = -0.25;
    out_dshape(2,1,0) =  0.25;
    out_dshape(3,1,0) = -0.25;

    out_dshape(0,1,1) = 0;
    out_dshape(1,1,1) = 0;
    out_dshape(2,1,1) = 0;
    out_dshape(3,1,1) = 0;
    break;

  case 8:
    out_dshape(0,0,0) = -0.5*(-1+eta);
    out_dshape(1,0,0) = -0.5*(-1+eta);
    out_dshape(2,0,0) =  0.5*( 1+eta);
    out_dshape(3,0,0) =  0.5*( 1+eta);
    out_dshape(4,0,0) = -1+eta;
    out_dshape(5,0,0) =  0.;
    out_dshape(6,0,0) = -1+eta;
    out_dshape(7,0,0) =  0.;

    out_dshape(0,0,1) = -0.25*(2*xi + eta) - 0.25*(-1+eta);
    out_dshape(1,0,1) =  0.25*(eta - 2*xi) + 0.25*(-1+eta);
    out_dshape(2,0,1) =  0.25*(2*xi + eta) + 0.25*( 1+eta);
    out_dshape(3,0,1) = -0.25*(eta - 2*xi) - 0.25*( 1+eta);
    out_dshape(4,0,1) =  xi;
    out_dshape(5,0,1) = -0.5*(-1+eta) - 0.5*(1+eta);
    out_dshape(6,0,1) = -xi;
    out_dshape(7,0,1) =  0.5*(-1+eta) + 0.5*(1+eta);

    out_dshape(0,1,0) = -0.25*(2*eta+xi) - 0.25*(-1+xi);
    out_dshape(1,1,0) =  0.25*(2*eta-xi) - 0.25*( 1+xi);
    out_dshape(2,1,0) =  0.25*(2*eta+xi) + 0.25*( 1+xi);
    out_dshape(3,1,0) = -0.25*(2*eta-xi) + 0.25*(-1+xi);
    out_dshape(4,1,0) =  0.5*(-1+xi) + 0.5*(1+xi);
    out_dshape(5,1,0) = -eta;
    out_dshape(6,1,0) = -0.5*(-1+xi) - 0.5*(1+xi);
    out_dshape(7,1,0) =  eta;

    out_dshape(0,1,1) = -0.5*(-1+xi);
    out_dshape(1,1,1) =  0.5*( 1+xi);
    out_dshape(2,1,1) =  0.5*( 1+xi);
    out_dshape(3,1,1) = -0.5*(-1+xi);
    out_dshape(4,1,1) =  0.;
    out_dshape(5,1,1) = -(1+xi);
    out_dshape(6,1,1) =  0.;
    out_dshape(7,1,1) =  (-1+xi);
    break;
  }
}

void shape_tri(const point &in_rs, vector<double> &out_shape)
{
  // NOTE: Reference triangle is defined in interval [0,1]^2

  // For the shape function for a general N-noded triangle, refer
  // to Finite Element Methods by Hughes, p. 166
  out_shape.resize(3); // nNodes

  out_shape[0] = in_rs.x;
  out_shape[1] = in_rs.y;
  out_shape[2] = 1 - in_rs.x - in_rs.y;
}

void shape_tri(const point &in_rs, double* out_shape)
{
  out_shape[0] = in_rs.x;
  out_shape[1] = in_rs.y;
  out_shape[2] = 1 - in_rs.x - in_rs.y;
}

void dshape_tri(point &, matrix<double> &out_dshape)
{
  out_dshape.setup(3,2); // nNodes, nDims

  out_dshape(0,0) =  1;
  out_dshape(1,0) =  0;
  out_dshape(2,0) = -1;

  out_dshape(0,1) =  0;
  out_dshape(1,1) =  1;
  out_dshape(2,1) = -1;
}

void shape_tet(const point &in_rs, vector<double> &out_shape)
{
  out_shape.resize(4); // nNodes
  shape_tet(in_rs, out_shape.data());
}

void shape_tet(const point &in_rs, double* out_shape)
{
  // NOTE: Reference tetrahedron is defined in interval [0,1]^3

  out_shape[0] = in_rs.x;
  out_shape[1] = in_rs.y;
  out_shape[2] = in_rs.z;
  out_shape[3] = 1 - in_rs.x - in_rs.y - in_rs.z;
}

void dshape_tet(point &, matrix<double> &out_dshape)
{
  out_dshape.setup(4,3); // nNodes, nDims

  out_dshape(0,0) =  1;
  out_dshape(1,0) =  0;
  out_dshape(2,0) =  0;
  out_dshape(3,0) = -1;

  out_dshape(0,1) =  0;
  out_dshape(1,1) =  1;
  out_dshape(2,1) =  0;
  out_dshape(3,1) = -1;

  out_dshape(0,2) =  0;
  out_dshape(1,2) =  0;
  out_dshape(2,2) =  1;
  out_dshape(3,2) = -1;
}

void getSimplex(int nDims, vector<double> x0, double L, matrix<double> X)
{
//  vector<double> dx(nDims+1,0);
//  X.setup(nDims,nDims+1);
//  X.initializeToZero();

//  for (int i=0; i<nDims; i++) {
//    dx[i+1] = sqrt(1-(i*dx[i]/(i+1))^2);
//    vector<int> ei(nDims); ei[i-1] = 1;
//    for (int j=0; j<nDims+1; j++) {
//      X(j,i+1) = 1/(i)*sum(X,2) + dx[i]*ei[j];
//    }
//  }

//  X *= L;
//  for (int i=0; i<nDims; i++) {
//    for (int j=0; j<nDims+1; j++) {
//      X(i,j) += x0[i];
//    }
//  }
}

vector<int> getOrder(vector<double> &data)
{
  vector<pair<double,size_t> > vp;
  vp.reserve(data.size());
  for (size_t i = 0 ; i != data.size() ; i++) {
    vp.push_back(make_pair(data[i], i));
  }

  // Sorting will put lower values [vp.first] ahead of larger ones,
  // resolving ties using the original index [vp.second]
  sort(vp.begin(), vp.end());
  vector<int> ind(data.size());
  for (size_t i = 0 ; i != vp.size() ; i++) {
    ind[i] = vp[i].second;
  }

  return ind;
}

Vec3 getFaceNormalTri(vector<point> &facePts, point &xc)
{
  point pt0 = facePts[0];
  point pt1 = facePts[1];
  point pt2 = facePts[2];
  Vec3 a = pt1 - pt0;
  Vec3 b = pt2 - pt0;
  Vec3 norm = a.cross(b);                         // Face normal vector
  Vec3 dxc = xc - (pt0+pt1+pt2)/3.;  // Face centroid to cell centroid
  if (norm*dxc > 0) {
    // Face normal is pointing into cell; flip
    norm /= -norm.norm();
  }
  else {
    // Face normal is pointing out of cell; keep direction
    norm /= norm.norm();
  }

  return norm;
}

Vec3 getFaceNormalQuad(vector<point> &facePts, point &xc)
{
  // Get the (approximate) face normal of an arbitrary 3D quadrilateral
  // by splitting into 2 triangles and averaging

  // Triangle #1
  point pt0 = facePts[0];
  point pt1 = facePts[1];
  point pt2 = facePts[2];
  Vec3 a = pt1 - pt0;
  Vec3 b = pt2 - pt0;
  Vec3 norm1 = a.cross(b);           // Face normal vector
  Vec3 dxc = xc - (pt0+pt1+pt2)/3.;  // Face centroid to cell centroid
  if (norm1*dxc > 0) {
    // Face normal is pointing into cell; flip
    norm1 *= -1;
  }

  // Triangle #2
  pt0 = facePts[3];
  a = pt1 - pt0;
  b = pt2 - pt0;
  Vec3 norm2 = a.cross(b);
  if (norm2*dxc > 0) {
    // Face normal is pointing into cell; flip
    norm2 *= -1;
  }

  // Average the two triangle's face outward normals
  Vec3 norm = (norm1+norm2)/2.;

  return norm;
}

Vec3 getEdgeNormal(vector<point> &edge, point &xc)
{
  Vec3 dx = edge[1] - edge[0];
  Vec3 norm = Vec3({-dx.y,dx.x,0});
  norm /= norm.norm();

  // Face centroid to cell centroid
  Vec3 dxc = xc - (edge[0]+edge[1])/2.;

  if (norm*dxc > 0) {
    // Face normal is pointing into cell; flip
    norm *= -1;
  }

  return norm;
}


void getBoundingBox(vector<point>& pts, point &minPt, point &maxPt)
{
  minPt = { INFINITY, INFINITY, INFINITY};
  maxPt = {-INFINITY,-INFINITY,-INFINITY};
  for (auto &pt:pts) {
    for (int dim=0; dim<3; dim++) {
      minPt[dim] = min(minPt[dim],pt[dim]);
      maxPt[dim] = max(maxPt[dim],pt[dim]);
    }
  }
}

void getBoundingBox(matrix<double>& pts, point &minPt, point &maxPt)
{
  minPt = { INFINITY, INFINITY, INFINITY};
  maxPt = {-INFINITY,-INFINITY,-INFINITY};
  for (int i=0; i<pts.getDim0(); i++) {
    for (int dim=0; dim<pts.getDim1(); dim++) {
      minPt[dim] = min(minPt[dim],pts(i,dim));
      maxPt[dim] = max(maxPt[dim],pts(i,dim));
    }
  }
}

void getBoundingBox(double *pts, int nPts, int nDims, point &minPt, point &maxPt)
{
  minPt = { INFINITY, INFINITY, INFINITY};
  maxPt = {-INFINITY,-INFINITY,-INFINITY};
  for (int i=0; i<nPts; i++) {
    for (int dim=0; dim<3; dim++) {
      minPt[dim] = min(minPt[dim],pts[i*nDims+dim]);
      maxPt[dim] = max(maxPt[dim],pts[i*nDims+dim]);
    }
  }
}

void getBoundingBox(double *pts, int nPts, int nDims, double *bbox)
{
  for (int i=0; i<nDims; i++) {
    bbox[i]       =  INFINITY;
    bbox[nDims+i] = -INFINITY;
  }

  for (int i=0; i<nPts; i++) {
    for (int dim=0; dim<nDims; dim++) {
      bbox[dim]       = min(bbox[dim],      pts[i*nDims+dim]);
      bbox[nDims+dim] = max(bbox[nDims+dim],pts[i*nDims+dim]);
    }
  }
}

vector<double> calcError(const vector<double> &U, const point &pos, input *params)
{
  if (params->testCase == 0) return U;

  int nDims = params->nDims;
  int nFields = params->nFields;

  vector<double> err(nFields);

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
      err[0] = rho;
      err[1] = rho * vx;
      err[2] = rho * vy;
      if (nDims == 3) err[3] = rho * vz;
      err[nDims+1] = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy + vz*vz));
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

      double x = fmod( (pos.x - params->time), (xmax - xmin) );
      double y = fmod( (pos.y - params->time), (ymax - ymin) );

      if (x > (xmax-xmin)/2) x -= (xmax-xmin);
      if (y > (ymax-ymin)/2) y -= (ymax-ymin);

      double f = 1.0 - (x*x + y*y);

      // Limiting rho to 1e-10 to avoid negative density/pressure issues
      double rho = max(pow(1. - eps*eps*(gamma-1.)/(8.*gamma*pi*pi)*exp(f), 1.0/(gamma-1.0) + 1e-5), 1e-10);
      double vx = 1. - eps*y / (2.*pi) * exp(f/2.);
      double vy = 1. + eps*x / (2.*pi) * exp(f/2.);
      double p = pow(rho,gamma);

      err[0] = rho;
      err[1] = rho * vx;
      err[2] = rho * vy;
      if (nDims == 3) err[3] = 0.;
      err[nDims+1] = p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
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

      double x = fmod( (pos.x - Uinf*cos(theta)*params->time), (xmax - xmin) );
      double y = fmod( (pos.y - Uinf*sin(theta)*params->time), (ymax - ymin) );

      if (x > (xmax-xmin)/2) x -= (xmax-xmin);
      if (y > (ymax-ymin)/2) y -= (ymax-ymin);

      double f = -(x*x + y*y) / (rc*rc);

      double vx = Uinf*(cos(theta) - y*eps/rc * exp(f/2.));
      double vy = Uinf*(sin(theta) + x*eps/rc * exp(f/2.));
      double rho = rhoInf*pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));
      double p   = Pinf  *pow(1. - (gamma-1.)/2. * eM * exp(f), gamma/(gamma-1.0));

      err[0] = rho;
      err[1] = rho * vx;
      err[2] = rho * vy;
      if (nDims == 3) err[3] = 0.;
      err[nDims+1] =  p/(gamma - 1) + (0.5*rho*(vx*vx + vy*vy));
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

    double x = abs(fmod( (pos.x - params->time), (xmax-xmin) ));
    double y = abs(fmod( (pos.y - params->time), (ymax-ymin) ));
    if (x > (xmax-xmin)/2) x -= (xmax-xmin);
    if (y > (ymax-ymin)/2) y -= (ymax-ymin);

    if (params->icType == 0) {
      /* --- Simple Gaussian bump centered at (0,0) --- */
      double r2 = x*x + y*y;
      err[0] = exp(-r2);
    }
    else if (params->icType == 1) {
      err[0] = 1 + sin(2.*pi*(pos.x+5.-params->time)/10.);
    }
    else if (params->icType == 2) {
      /* --- Test case for debugging - cos(x)*cos(y)*cos(z) over domain --- */
      err[0] = cos(2*pi*pos.x/6.)*cos(2*pi*pos.y/6.)*cos(2*pi*pos.z/6.);
    }
  }

  for (int i=0; i<nFields; i++)
    err[i] = U[i] - err[i];

  if (params->errorNorm == 1)
    for (auto &val:err) val = abs(val); // L1 norm
  else if (params->errorNorm == 2)
    for (auto &val:err) val *= val;     // L2 norm

  return err;
}
