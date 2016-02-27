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

#include <cmath>
#include <set>

#include "flux.hpp"
#include "geo.hpp"
#include "global.hpp"
#include "polynomials.hpp"

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
  if (nNodes == 4) {
      out_shape[0] = 0.25*(1-xi)*(1-eta);
      out_shape[1] = 0.25*(1+xi)*(1-eta);
      out_shape[2] = 0.25*(1+xi)*(1+eta);
      out_shape[3] = 0.25*(1-xi)*(1+eta);
  }
  else if (nNodes == 8) {
      out_shape[0] = -0.25*(1-xi)*(1-eta)*(1+eta+xi);
      out_shape[1] = -0.25*(1+xi)*(1-eta)*(1+eta-xi);
      out_shape[2] = -0.25*(1+xi)*(1+eta)*(1-eta-xi);
      out_shape[3] = -0.25*(1-xi)*(1+eta)*(1-eta+xi);
      out_shape[4] = 0.5*(1-xi)*(1+ xi)*(1-eta);
      out_shape[5] = 0.5*(1+xi)*(1+eta)*(1-eta);
      out_shape[6] = 0.5*(1-xi)*(1+ xi)*(1+eta);
      out_shape[7] = 0.5*(1-xi)*(1+eta)*(1-eta);
  }
  else {
    int nSide = sqrt(nNodes);

    if (nSide*nSide != nNodes)
      FatalError("For Lagrange quad of order N, must have (N+1)^2 shape points.");

    vector<double> xlist(nSide);
    double dxi = 2./(nSide-1);

    for (int i=0; i<nSide; i++)
      xlist[i] = -1. + i*dxi;

    int nLevels = nSide / 2;
    int isOdd = nSide % 2;

    /* Recursion for all high-order Lagrange elements:
     * 4 corners, each edge's points, interior points */
    int nPts = 0;
    for (int i = 0; i < nLevels; i++) {
      // Corners
      int i2 = (nSide-1) - i;
      out_shape[nPts+0] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i);
      out_shape[nPts+1] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i);
      out_shape[nPts+2] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2);
      out_shape[nPts+3] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2);
      nPts += 4;

      // Edges: Bottom, right, top, left
      int nSide2 = nSide - 2 * (i+1);
      for (int j = 0; j < nSide2; j++) {
        out_shape[nPts+0*nSide2+j] = Lagrange(xlist, xi, i+1+j) * Lagrange(xlist, eta, i);
        out_shape[nPts+1*nSide2+j] = Lagrange(xlist, xi, i2) * Lagrange(xlist, eta, i+1+j);
        out_shape[nPts+2*nSide2+j] = Lagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2);
        out_shape[nPts+3*nSide2+j] = Lagrange(xlist, xi, i) * Lagrange(xlist, eta, i2-1-j);
      }
      nPts += 4*nSide2;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd) {
      out_shape[nNodes-1] = Lagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2);
    }
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
  out_dshape.setup(nNodes,2);
  dshape_quad(in_rs, &out_dshape(0,0), nNodes);
}

void dshape_quad(const point &in_rs, double* out_dshape, int nNodes)
{
  double xi  = in_rs.x;
  double eta = in_rs.y;

  if (nNodes == 4) {
    out_dshape[0*2+0] = -0.25*(1-eta);
    out_dshape[1*2+0] =  0.25*(1-eta);
    out_dshape[2*2+0] =  0.25*(1+eta);
    out_dshape[3*2+0] = -0.25*(1+eta);

    out_dshape[0*2+1] = -0.25*(1-xi);
    out_dshape[1*2+1] = -0.25*(1+xi);
    out_dshape[2*2+1] =  0.25*(1+xi);
    out_dshape[3*2+1] =  0.25*(1-xi);
  }
  else if (nNodes == 8) {
    out_dshape[0*2+0] = -0.25*(-1+eta)*(2*xi+eta);
    out_dshape[1*2+0] =  0.25*(-1+eta)*(eta - 2*xi);
    out_dshape[2*2+0] =  0.25*( 1+eta)*(2*xi+eta);
    out_dshape[3*2+0] = -0.25*( 1+eta)*(eta-2*xi);
    out_dshape[4*2+0] =    xi*(-1+eta);
    out_dshape[5*2+0] = -0.5 *( 1+eta)*(-1+eta);
    out_dshape[6*2+0] =   -xi*( 1+eta);
    out_dshape[7*2+0] =  0.5 *( 1+eta)*(-1+eta);

    out_dshape[0*2+1] = -0.25*(-1+xi)*(2*eta+xi);
    out_dshape[1*2+1] =  0.25*( 1+xi)*(2*eta - xi);
    out_dshape[2*2+1] =  0.25*( 1+xi)*(2*eta+xi);
    out_dshape[3*2+1] = -0.25*(-1+xi)*(2*eta-xi);
    out_dshape[4*2+1] =  0.5 *( 1+xi)*(-1+xi);
    out_dshape[5*2+1] =  -eta*( 1+xi);
    out_dshape[6*2+1] = -0.5 *( 1+xi)*(-1+xi);
    out_dshape[7*2+1] =   eta*(-1+xi);
  }
  else {
    int nSide = sqrt(nNodes);

    if (nSide*nSide != nNodes)
      FatalError("For Lagrange quad of order N, must have (N+1)^2 shape points.");

    vector<double> xlist(nSide);
    double dxi = 2./(nSide-1);

    for (int i=0; i<nSide; i++)
      xlist[i] = -1. + i*dxi;

    int nLevels = nSide / 2;
    int isOdd = nSide % 2;

    /* Recursion for all high-order Lagrange elements:
     * 4 corners, each edge's points, interior points */
    int nPts = 0;
    for (int i = 0; i < nLevels; i++) {
      // Corners
      int i2 = (nSide-1) - i;
      out_dshape[nPts+0] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i);
      out_dshape[nPts+1] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i);
      out_dshape[nPts+2] = dLagrange(xlist, xi, i2) * Lagrange(xlist, eta, i2);
      out_dshape[nPts+3] = dLagrange(xlist, xi, i)  * Lagrange(xlist, eta, i2);

      out_dshape[2*(nPts+0)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i);
      out_dshape[2*(nPts+1)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i);
      out_dshape[2*(nPts+2)+1] = Lagrange(xlist, xi, i2) * dLagrange(xlist, eta, i2);
      out_dshape[2*(nPts+3)+1] = Lagrange(xlist, xi, i)  * dLagrange(xlist, eta, i2);
      nPts += 4;

      // Edges
      int nSide2 = nSide - 2 * (i+1);
      for (int j = 0; j < nSide2; j++) {
        out_dshape[2*(nPts+0*nSide2+j)] = dLagrange(xlist, xi, i+1+j)  * Lagrange(xlist, eta, i);
        out_dshape[2*(nPts+1*nSide2+j)] = dLagrange(xlist, xi, i2)   * Lagrange(xlist, eta, i+1+j);
        out_dshape[2*(nPts+2*nSide2+j)] = dLagrange(xlist, xi, i2-1-j) * Lagrange(xlist, eta, i2);
        out_dshape[2*(nPts+3*nSide2+j)] = dLagrange(xlist, xi, i)    * Lagrange(xlist, eta, i2-1-j);

        out_dshape[2*(nPts+0*nSide2+j)+1] = Lagrange(xlist, xi, i+1+j)  * dLagrange(xlist, eta, i);
        out_dshape[2*(nPts+1*nSide2+j)+1] = Lagrange(xlist, xi, i2)   * dLagrange(xlist, eta, i+1+j);
        out_dshape[2*(nPts+2*nSide2+j)+1] = Lagrange(xlist, xi, i2-1-j) * dLagrange(xlist, eta, i2);
        out_dshape[2*(nPts+3*nSide2+j)+1] = Lagrange(xlist, xi, i)    * dLagrange(xlist, eta, i2-1-j);
      }
      nPts += 4*nSide2;
    }

    // Center node for even-ordered Lagrange quads (odd value of nSide)
    if (isOdd) {
      out_dshape[2*(nNodes-1)+0] = dLagrange(xlist, xi, nSide/2) * Lagrange(xlist, eta, nSide/2);
      out_dshape[2*(nNodes-1)+1] = Lagrange(xlist, xi, nSide/2) * dLagrange(xlist, eta, nSide/2);
    }
  }
}

void dshape_hex(const point &in_rst, matrix<double> &out_dshape, int nNodes)
{
  out_dshape.setup(nNodes,3);
  dshape_hex(in_rst, &out_dshape(0,0), nNodes);
}

void dshape_hex(const point &in_rst, double* out_dshape, int nNodes)
{
  double xi  = in_rst.x;
  double eta = in_rst.y;
  double mu = in_rst.z;
  switch(nNodes) {
    case 8:
      out_dshape[3*0+0] = -0.125*(1-eta)*(1-mu);
      out_dshape[3*1+0] =  0.125*(1-eta)*(1-mu);
      out_dshape[3*2+0] =  0.125*(1+eta)*(1-mu);
      out_dshape[3*3+0] = -0.125*(1+eta)*(1-mu);

      out_dshape[3*4+0] = -0.125*(1-eta)*(1+mu);
      out_dshape[3*5+0] =  0.125*(1-eta)*(1+mu);
      out_dshape[3*6+0] =  0.125*(1+eta)*(1+mu);
      out_dshape[3*7+0] = -0.125*(1+eta)*(1+mu);

      out_dshape[3*0+1] = -0.125*(1-xi)*(1-mu);
      out_dshape[3*1+1] = -0.125*(1+xi)*(1-mu);
      out_dshape[3*2+1] =  0.125*(1+xi)*(1-mu);
      out_dshape[3*3+1] =  0.125*(1-xi)*(1-mu);

      out_dshape[3*4+1] = -0.125*(1-xi)*(1+mu);
      out_dshape[3*5+1] = -0.125*(1+xi)*(1+mu);
      out_dshape[3*6+1] =  0.125*(1+xi)*(1+mu);
      out_dshape[3*7+1] =  0.125*(1-xi)*(1+mu);

      out_dshape[3*0+2] = -0.125*(1-xi)*(1-eta);
      out_dshape[3*1+2] = -0.125*(1+xi)*(1-eta);
      out_dshape[3*2+2] = -0.125*(1+xi)*(1+eta);
      out_dshape[3*3+2] = -0.125*(1-xi)*(1+eta);

      out_dshape[3*4+2] =  0.125*(1-xi)*(1-eta);
      out_dshape[3*5+2] =  0.125*(1+xi)*(1-eta);
      out_dshape[3*6+2] =  0.125*(1+xi)*(1+eta);
      out_dshape[3*7+2] =  0.125*(1-xi)*(1+eta);
      break;
    case 20: {
      double XI[8]  = {-1,1,1,-1,-1,1,1,-1};
      double ETA[8] = {-1,-1,1,1,-1,-1,1,1};
      double MU[8]  = {-1,-1,-1,-1,1,1,1,1};
      // Corner Nodes
      for (int i=0; i<8; i++) {
        out_dshape[3*i+0] = .125*XI[i] *(1+eta*ETA[i])*(1 + mu*MU[i])*(2*xi*XI[i] +   eta*ETA[i] +   mu*MU[i]-1);
        out_dshape[3*i+1] = .125*ETA[i]*(1 + xi*XI[i])*(1 + mu*MU[i])*(  xi*XI[i] + 2*eta*ETA[i] +   mu*MU[i]-1);
        out_dshape[3*i+2] = .125*MU[i] *(1 + xi*XI[i])*(1+eta*ETA[i])*(  xi*XI[i] +   eta*ETA[i] + 2*mu*MU[i]-1);
      }
      // Edge Nodes, xi = 0
      out_dshape[ 3*8+0] = -.5*xi*(1-eta)*(1-mu);  out_dshape[ 3*8+1] = -.25*(1-xi*xi)*(1-mu);  out_dshape[ 3*8+2] = -.25*(1-xi*xi)*(1-eta);
      out_dshape[3*10+0] = -.5*xi*(1+eta)*(1-mu);  out_dshape[3*10+1] =  .25*(1-xi*xi)*(1-mu);  out_dshape[3*10+2] = -.25*(1-xi*xi)*(1+eta);
      out_dshape[3*16+0] = -.5*xi*(1-eta)*(1+mu);  out_dshape[3*16+1] = -.25*(1-xi*xi)*(1+mu);  out_dshape[3*16+2] =  .25*(1-xi*xi)*(1-eta);
      out_dshape[3*18+0] = -.5*xi*(1+eta)*(1+mu);  out_dshape[3*18+1] =  .25*(1-xi*xi)*(1+mu);  out_dshape[3*18+2] =  .25*(1-xi*xi)*(1+eta);
      // Edge Nodes, eta = 0
      out_dshape[ 3*9+1] = -.5*eta*(1+xi)*(1-mu);  out_dshape[ 3*9+0] =  .25*(1-eta*eta)*(1-mu);  out_dshape[ 3*9+2] = -.25*(1-eta*eta)*(1+xi);
      out_dshape[3*11+1] = -.5*eta*(1-xi)*(1-mu);  out_dshape[3*11+0] = -.25*(1-eta*eta)*(1-mu);  out_dshape[3*11+2] = -.25*(1-eta*eta)*(1-xi);
      out_dshape[3*17+1] = -.5*eta*(1+xi)*(1+mu);  out_dshape[3*17+0] =  .25*(1-eta*eta)*(1+mu);  out_dshape[3*17+2] =  .25*(1-eta*eta)*(1+xi);
      out_dshape[3*19+1] = -.5*eta*(1-xi)*(1+mu);  out_dshape[3*19+0] = -.25*(1-eta*eta)*(1+mu);  out_dshape[3*19+2] =  .25*(1-eta*eta)*(1-xi);
      // Edge Nodes, mu = 0;
      out_dshape[3*12+2] = -.5*mu*(1-xi)*(1-eta);  out_dshape[3*12+0] = -.25*(1-mu*mu)*(1-eta);  out_dshape[3*12+1] = -.25*(1-mu*mu)*(1-xi);
      out_dshape[3*13+2] = -.5*mu*(1+xi)*(1-eta);  out_dshape[3*13+0] =  .25*(1-mu*mu)*(1-eta);  out_dshape[3*13+1] = -.25*(1-mu*mu)*(1+xi);
      out_dshape[3*14+2] = -.5*mu*(1+xi)*(1+eta);  out_dshape[3*14+0] =  .25*(1-mu*mu)*(1+eta);  out_dshape[3*14+1] =  .25*(1-mu*mu)*(1+xi);
      out_dshape[3*15+2] = -.5*mu*(1-xi)*(1+eta);  out_dshape[3*15+0] = -.25*(1-mu*mu)*(1+eta);  out_dshape[3*15+1] =  .25*(1-mu*mu)*(1-xi);
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

vector<double> calcError(const double* const U, const point &pos, input *params)
{
  if (params->testCase == 0) {
    vector<double> U_out(U,U+params->nFields);
    return U_out;
  }

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
      if (x > xmax) x -= (xmax-xmin);
      if (y > ymax) y -= (ymax-ymin);
      if (x < xmin) x += (xmax-xmin);
      if (y < ymin) y += (ymax-ymin);

      double f = 1.0 - (x*x + y*y);

      // Limiting rho to 1e-10 to avoid negative density/pressure issues
      double rho = max(pow(1. - eps*eps*(gamma-1.)/(8.*gamma*pi*pi)*exp(f), 1.0/(gamma-1.0)), 1e-10);
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

      double advX = Uinf*cos(theta)*params->time;
      double advY = Uinf*sin(theta)*params->time;
      double x = (fmod( (pos.x - advX), (xmax - xmin) ));
      double y = (fmod( (pos.y - advY), (ymax - ymin) ));
      if (x > xmax) x -= (xmax-xmin);
      if (x < xmin) x += (xmax-xmin);
      if (y > ymax) y -= (ymax-ymin);
      if (y < ymin) y += (ymax-ymin);

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
    for (auto &val:err) val = std::abs(val); // L1 norm
  else if (params->errorNorm == 2)
    for (auto &val:err) val *= val;     // L2 norm

  return err;
}

void calcSolutionFromFlux(matrix<double> &F, vector<double> &U, input *params)
{
  // Use Newton iterations to converge to U(F) such that F(U) = F
  int nFields = params->nFields;
  int nDims = params->nDims;

  matrix<double> A,B;
  matrix<double> tempF(nDims,nFields);
  vector<double> dF(nFields), dG(nFields), dH;
  if (nDims==3) dH.resize(nFields);

  inviscidFlux(U.data(),tempF,params);

  for (int i=0; i<nFields; i++) {
    dF[i] = F(0,i) - tempF(0,i);
    dG[i] = F(1,i) - tempF(1,i);
    if (nDims==3) dH[i] = F(2,i) - tempF(2,i);
  }

  // Magnitude of difference between current and desired flux vector
  double norm = 0.;
  for (auto &val:dF) norm += val*val;
  for (auto &val:dG) norm += val*val;
  for (auto &val:dH) norm += val*val;

  // Hackish, yeah...
  vector<double> maxDU = U;
  for (auto &val:maxDU) val = std::max(.25*std::abs(val),.1);

  int iter = 0;
  double tol = 1e-10;
  if (nDims == 2) {
    while (norm > tol) {
      calcFluxJacobian2D(U,A,B,params);

      auto du1 = A.solve(dF);
      auto du2 = B.solve(dG);

      for (int i=0; i<nFields; i++) {
        U[i] += (du1[i]+du2[i])/2.;
      }

      inviscidFlux(U.data(),tempF,params);

      for (int i=0; i<nFields; i++) {
        dF[i] = F(0,i) - tempF(0,i);
        dG[i] = F(1,i) - tempF(1,i);
      }

      norm = 0.;
      for (auto &val:dF) norm += val*val;
      for (auto &val:dG) norm += val*val;

      iter++;
      if (iter>10) break;
    }
  }
  else {

  }
}

void calcFluxJacobian2D(const vector<double> &U, matrix<double> &dFdU, matrix<double> &dGdU, input *params)
{
  int nFields = U.size();
  if (nFields != 4) FatalError("flux Jacobian only implemented for 2D Euler.");

  dFdU.setup(4,4);
  dGdU.setup(4,4);

  double gamma = params->gamma;

  double rho = U[0];
  double u = U[1]/rho;
  double v = U[2]/rho;
  double uv2 = u*u+v*v;
  double e = U[3];
  double E = e/rho;
  double p = (gamma-1.0)*(U[3]-rho*(u*u+v*v));

  double a0 = (gamma-1)*(u*u+v*v);
  double a1 = a0 - gamma*e/rho;

  // A = dFdU
  dFdU(0,1) = 1;

  dFdU(1,0) = a0/2. - u*u;
  dFdU(1,1) = (3-gamma)*u;
  dFdU(1,2) = (1-gamma)*v;
  dFdU(1,3) = gamma-1;

  dFdU(2,0) = -u*v;
  dFdU(2,1) = v;
  dFdU(2,2) = u;

  dFdU(3,0) = a1*u;
  dFdU(3,1) = gamma*e/rho - a0/2. - (gamma-1)*u*u;
  dFdU(3,2) = (1-gamma)*u*v;
  dFdU(3,3) = gamma*u;

  // B = dGdU
  dGdU(0,2) = 1;

  dGdU(1,0) = -u*v;
  dGdU(1,1) = v;
  dGdU(1,2) = u;

  dGdU(2,0) = a0/2. - v*v;
  dGdU(2,1) = (1-gamma)*u;
  dGdU(2,2) = (3-gamma)*v;
  dGdU(2,3) = gamma-1;

  dGdU(3,0) = a1*v;
  dGdU(3,1) = (1-gamma)*u*v;
  dGdU(3,2) = gamma*e/rho - a0/2. - (gamma-1)*v*v;
  dGdU(3,3) = gamma*v;
}

void refineGrid2D(geo &grid_c, geo &grid_f, int nLevels, int nNodes_c, int shapeOrder_f)
{
  int nCellSplit = pow(4, nLevels);
  int nFaceSplit = pow(2, nLevels);

  int shapeOrder_c;
  if (nNodes_c == 4)
    shapeOrder_c = 1;
  else if (nNodes_c == 8)
    shapeOrder_c = 2;
  else
    shapeOrder_c = std::sqrt(nNodes_c) - 1;

  if (shapeOrder_c < 5 && grid_c.rank == 0)
    cout << "WARNING: For best results when refining grids, use a very high order shape representation on the coarse grid!" << endl;

  int nSide_c = shapeOrder_c + 1;
  int nSide_f = shapeOrder_f + 1;
  int nNodes_f = nSide_f * nSide_f;

  int nEles_c = grid_c.c2v.getDim0();
  int nEles_f = nEles_c * nCellSplit;

  int nFaces_c = grid_c.f2v.getDim0();
  int nFaces_f = nFaces_c * nFaceSplit
               + nEles_c * 4 * (nCellSplit - 1.) / (3.);

  grid_f.nEles = nEles_f;
  grid_f.c2v.setup(nEles_f, nNodes_f);

  vector<point> xv_new;

  /* Setup arrays for introducing new nodes inside each coarse-grid cell */
  int ndSplit1D = shapeOrder_f * nFaceSplit + 1;
  double dxi_nd = 2. / (ndSplit1D - 1.);
  vector<double> xiList(ndSplit1D);
  for (int i = 0; i < ndSplit1D; i++)
    xiList[i] = -1. + i * dxi_nd;

  for (int ic = 0; ic < nEles_c; ic++) {
    int nv_cur = xv_new.size();

    /* Get physical position of new nodes */
    vector<double> shape_tmp;
    for (int j = 0; j < ndSplit1D; j++) {
      for (int i = 0; i < ndSplit1D; i++) {
        point loc = {xiList[i], xiList[j], 0.};
        point pos;
        shape_quad(loc, shape_tmp, nNodes_c);
        for (int n = 0; n < nNodes_c; n++) {
          pos.x += shape_tmp[n] * grid_c.xv(grid_c.c2v(ic,n), 0);
          pos.y += shape_tmp[n] * grid_c.xv(grid_c.c2v(ic,n), 1);
        }
        xv_new.push_back(pos);
      }
    }

    /* Setup connectivity of new fine-grid sub-cells */
    for (int i = 0; i < nFaceSplit; i++) {
      for (int j = 0; j < nFaceSplit; j++) {
        // i,j coords of 'bottom-left' node in new cell
        int ioff = i * shapeOrder_f;
        int joff = j * shapeOrder_f;
        int ic_new = ic*nCellSplit + i * (nFaceSplit) + j;

        /* Recursion for high-order shape functions:
         * 4 corners, each edge's points, interior points
         * i is 'row', j is 'col', starting from bottom-left */
        int nSLevels = nSide_f / 2;
        int isOdd = nSide_f % 2;

        int nPts = 0;
        int start = nv_cur + ioff * ndSplit1D + joff;

        auto getCoord = [=](int row, int col) -> int {
          return start + row * ndSplit1D + col;
        };

        for (int ii = 0; ii < nSLevels; ii++) {
          // Corners
          int i2 = (nSide_f-1) - ii;
          grid_f.c2v(ic_new, nPts+0) = getCoord(ii, ii);
          grid_f.c2v(ic_new, nPts+1) = getCoord(ii, i2);
          grid_f.c2v(ic_new, nPts+2) = getCoord(i2, i2);
          grid_f.c2v(ic_new, nPts+3) = getCoord(i2, ii);
          nPts += 4;

          // Edges: Bottom, right, top, left
          int nSide2 = nSide_f - 2 * (ii+1);
          for (int jj = 0; jj < nSide2; jj++) {
            grid_f.c2v(ic_new, nPts+0*nSide2+jj) = getCoord(ii, ii+1+jj);
            grid_f.c2v(ic_new, nPts+1*nSide2+jj) = getCoord(ii+1+jj, i2);
            grid_f.c2v(ic_new, nPts+2*nSide2+jj) = getCoord(i2, i2-1-jj);
            grid_f.c2v(ic_new, nPts+3*nSide2+jj) = getCoord(i2-1-jj, ii);
          }
          nPts += 4*nSide2;
        }

        // Center node for even-ordered Lagrange quads (odd value of nSide)
        if (isOdd) {
          grid_f.c2v(ic_new, nNodes_f-1) = getCoord(nSide_f/2, nSide_f/2);
        }
      }
    }
  }


  /* Use ptMap to redirect higher-numbered duplicate node to lower-numbered node */

  int nVerts_f = xv_new.size();
  vector<vector<int>> boundPoints(grid_c.nBounds);
  vector<int> ptMap(nVerts_f);
  for (int iv = 0; iv < nVerts_f; iv++) ptMap[iv] = iv;

  grid_f.bcList = grid_c.bcList;
  grid_f.nBounds = grid_c.nBounds;
  grid_f.nBndPts.resize(grid_f.nBounds);

  for (int F = 0; F < nFaces_c; F++) {
    int ic1 = grid_c.f2c(F,0);
    int ic2 = grid_c.f2c(F,1);
    if (ic2 >= 0) {
      // Screw it - just do brute-force comparison
      int start1 = ic1 * ndSplit1D * ndSplit1D;
      int start2 = ic2 * ndSplit1D * ndSplit1D;
      for (int i = 0; i < ndSplit1D*ndSplit1D; i++) {
        int iv1 = start1 + i;
        point pt1 = xv_new[iv1];
        for (int j = 0; j < ndSplit1D*ndSplit1D; j++) {
          int iv2 = start2 + j;
          Vec3 D = xv_new[iv2] - pt1;
          double dist = D.norm();
          if (std::abs(dist) < 1e-12) {
            if (iv1 > iv2)
              ptMap[iv1] = ptMap[iv2];
            else
              ptMap[iv2] = ptMap[iv1];
          }
        }
      }
    }
    else {
      // Boundary face - figure out which bnd and add points to BC pt list
      auto cellFaces = grid_c.c2f.getRow(ic1);
      int locF = findFirst(cellFaces,F);
      int bfID = findFirst(grid_c.bndFaces,F);
      if (bfID >= 0) {
        int BC = grid_c.bcType[bfID];
        int bcid = findFirst(grid_c.bcList, BC);

        int start, stride;
        int ic_start = ic1 * ndSplit1D * ndSplit1D;
        if (locF == 0) {
          // Bottom
          start = ic_start;
          stride = 1;
        }
        else if (locF == 1) {
          // Right
          start = ic_start + ndSplit1D - 1;
          stride = ndSplit1D;
        }
        else if (locF == 2) {
          // Top
          start = ic_start + ndSplit1D * ndSplit1D - 1;
          stride = -1;
        }
        else if (locF == 3) {
          // Left
          start = ic_start + ndSplit1D * (ndSplit1D - 1);
          stride = -ndSplit1D;
        }
        else
          FatalError("Improper locF value.");

        for (int i = 0; i < ndSplit1D; i++) {
          boundPoints[bcid].push_back(start + i * stride);
        }
      }
    }
  }

  /* Remove duplicated nodes */

  for (int iv = nVerts_f - 1; iv >= 0; iv--) {
    if (ptMap[iv] < iv) {
      xv_new.erase(xv_new.begin()+iv);
      for (int iv2 = iv + 1; iv2 < nVerts_f; iv2++) {
        if (ptMap[iv2] > iv)
          ptMap[iv2]--;
      }
    }
  }

  /* Copy nodes into grid */

  grid_f.nVerts = xv_new.size();
  grid_f.xv.setup(xv_new.size(), 2);
  for (int iv = 0; iv < grid_f.nVerts; iv++) {
    grid_f.xv(iv,0) = xv_new[iv].x;
    grid_f.xv(iv,1) = xv_new[iv].y;
  }

  /* Update connectivity */

  for (int ic = 0; ic < nEles_f; ic++)
    for (int j = 0; j < nNodes_f; j++)
      grid_f.c2v(ic,j) = ptMap[grid_f.c2v(ic,j)];

  /* Setup boundary-condition data */

  for (auto &bndVec:boundPoints) {
    for (auto &iv:bndVec)
      iv = ptMap[iv];

    std::sort(bndVec.begin(), bndVec.end());
    auto it = std::unique(bndVec.begin(), bndVec.end());
    bndVec.resize(std::distance(bndVec.begin(), it));
  }

  for (int bc = 0; bc < grid_f.nBounds; bc++)
    grid_f.nBndPts[bc] = boundPoints[bc].size();

  int maxNBndPts = getMax(grid_f.nBndPts);
  grid_f.bndPts.setup(grid_f.nBounds,maxNBndPts);
  for (int bc = 0; bc < grid_f.nBounds; bc++) {
    for (int j = 0; j < boundPoints[bc].size(); j++) {
      grid_f.bndPts(bc,j) = boundPoints[bc][j];
    }
  }

  /* Other */

  grid_f.ctype.assign(nEles_f, QUAD);
  grid_f.c2nv.assign(nEles_f, nNodes_f);
  grid_f.c2nf.assign(nEles_f, 4);
}
