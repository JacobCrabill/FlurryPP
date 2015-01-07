/*!
 * \file polynomials.cpp
 * \brief Polynomial definitions for all Flux Reconstruction operations
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/polynomials.hpp"

#include <cmath>

using namespace std;


double Lagrange(vector<double> &x_lag, double &y, unsigned int &mode)
{
  unsigned int i;
  double lag = 1.0;

  for(i=0; i<x_lag.size(); i++) {
    if(i!=mode) {
      lag = lag*((y-x_lag[i])/(x_lag[mode]-x_lag[i]));
    }
  }

  return lag;
}


double dLagrange(vector<double> &x_lag, double &y, unsigned int &mode)
{
  unsigned int i, j;
  double dLag, dLag_num, dLag_den;

  dLag = 0.0;

  for (i=0; i<x_lag.size(); i++) {
    if (i!=mode) {
      dLag_num = 1.0;
      dLag_den = 1.0;

      for (j=0; j<x_lag.size(); j++) {
        if (j!=mode && j!=i) {
          dLag_num = dLag_num*(y-x_lag[j]);
        }

        if (j!=mode) {
          dLag_den = dLag_den*(x_lag[mode]-x_lag[j]);
        }
      }

      dLag = dLag+(dLag_num/dLag_den);
    }
  }

  return dLag;
}

double ddLagrange(vector<double> &x_lag, double &y, unsigned int &mode)
{
  unsigned int i, j, k;
  double ddLag, ddLag_num, ddLag_den;

  ddLag = 0.0;

  for (i=0; i<x_lag.size(); i++) {
    if (i!=mode) {

      for (j=0; j<x_lag.size(); j++) {
        if (j!=mode) {
          ddLag_num = 1.0;
          ddLag_den = 1.0;

          for (k=0; k<x_lag.size(); k++) {
            if(k!=mode) {
              if(k!=i && k!=j) {
                ddLag_num *= (y-x_lag[k]);
              }

              ddLag_den *= (x_lag[mode]-x_lag[k]);
            }
          }

          if(j!=i) {
            ddLag = ddLag + (ddLag_num/ddLag_den);
          }
        }
      }

    }
  }

  return ddLag;
}


void shape_quad(point &in_rs, vector<double> &out_shape)
{
  out_shape.resize(4); // nNodes

  out_shape[0] = 0.25*(1-in_rs.x)*(1-in_rs.y);
  out_shape[1] = 0.25*(1+in_rs.x)*(1-in_rs.y);
  out_shape[2] = 0.25*(1+in_rs.x)*(1+in_rs.y);
  out_shape[3] = 0.25*(1-in_rs.x)*(1+in_rs.y);
}

void dshape_quad(point &in_rs, matrix<double> &out_dshape)
{
  out_dshape.setup(4,2);

  out_dshape[0][0] = -0.25*(1-in_rs.y);
  out_dshape[1][0] =  0.25*(1-in_rs.y);
  out_dshape[2][0] =  0.25*(1+in_rs.y);
  out_dshape[3][0] = -0.25*(1+in_rs.y);

  out_dshape[0][1] = -0.25*(1-in_rs.x);
  out_dshape[1][1] = -0.25*(1+in_rs.x);
  out_dshape[2][1] =  0.25*(1+in_rs.x);
  out_dshape[3][1] =  0.25*(1-in_rs.x);
}

void shape_tri(point &in_rs, vector<double> &out_shape)
{
  // For the shape function for a general N-noded triangle, refer
  // to Finite Element Methods by Hughes, p. 166
  out_shape.resize(3); // nNodes

  out_shape[0] = in_rs.x;
  out_shape[1] = in_rs.y;
  out_shape[2] = 1 - in_rs.x - in_rs.y;
}

void dshape_tri(point &in_rs, matrix<double> &out_dshape)
{
  out_dshape.setup(3,2); // nNodes, nDims

  out_dshape[0][0] =  1;
  out_dshape[1][0] =  0;
  out_dshape[2][0] = -1;

  out_dshape[0][1] =  0;
  out_dshape[1][1] =  1;
  out_dshape[2][1] = -1;
}

// helper method to evaluate a normalized jacobi polynomial
double eval_jacobi(double in_r, int in_alpha, int in_beta, int in_mode)
{
  double jacobi;

  if(in_mode==0) {
    double dtemp_0, dtemp_1, dtemp_2;

    dtemp_0=pow(2.0,(-in_alpha-in_beta-1));
    dtemp_1=eval_gamma(in_alpha+in_beta+2);
    dtemp_2=eval_gamma(in_alpha+1)*eval_gamma(in_beta+1);

    jacobi=sqrt(dtemp_0*(dtemp_1/dtemp_2));
  }
  else if (in_mode==1) {
    double dtemp_0, dtemp_1, dtemp_2, dtemp_3, dtemp_4, dtemp_5;

    dtemp_0=pow(2.0,(-in_alpha-in_beta-1));
    dtemp_1=eval_gamma(in_alpha+in_beta+2);
    dtemp_2=eval_gamma(in_alpha+1)*eval_gamma(in_beta+1);
    dtemp_3=in_alpha+in_beta+3;
    dtemp_4=(in_alpha+1)*(in_beta+1);
    dtemp_5=(in_r*(in_alpha+in_beta+2)+(in_alpha-in_beta));

    jacobi=0.5*sqrt(dtemp_0*(dtemp_1/dtemp_2))*sqrt(dtemp_3/dtemp_4)*dtemp_5;
  }
  else {
    double dtemp_0, dtemp_1, dtemp_3, dtemp_4, dtemp_5, dtemp_6, dtemp_7, dtemp_8, dtemp_9, dtemp_10, dtemp_11, dtemp_12, dtemp_13, dtemp_14;

    dtemp_0=in_mode*(in_mode+in_alpha+in_beta)*(in_mode+in_alpha)*(in_mode+in_beta);
    dtemp_1=((2*in_mode)+in_alpha+in_beta-1)*((2*in_mode)+in_alpha+in_beta+1);
    dtemp_3=(2*in_mode)+in_alpha+in_beta;

    dtemp_4=(in_mode-1)*((in_mode-1)+in_alpha+in_beta)*((in_mode-1)+in_alpha)*((in_mode-1)+in_beta);
    dtemp_5=((2*(in_mode-1))+in_alpha+in_beta-1)*((2*(in_mode-1))+in_alpha+in_beta+1);
    dtemp_6=(2*(in_mode-1))+in_alpha+in_beta;

    dtemp_7=-((in_alpha*in_alpha)-(in_beta*in_beta));
    dtemp_8=((2*(in_mode-1))+in_alpha+in_beta)*((2*(in_mode-1))+in_alpha+in_beta+2);

    dtemp_9=(2.0/dtemp_3)*sqrt(dtemp_0/dtemp_1);
    dtemp_10=(2.0/dtemp_6)*sqrt(dtemp_4/dtemp_5);
    dtemp_11=dtemp_7/dtemp_8;

    dtemp_12=in_r*eval_jacobi(in_r,in_alpha,in_beta,in_mode-1);
    dtemp_13=dtemp_10*eval_jacobi(in_r,in_alpha,in_beta,in_mode-2);
    dtemp_14=dtemp_11*eval_jacobi(in_r,in_alpha,in_beta,in_mode-1);

    jacobi=(1.0/dtemp_9)*(dtemp_12-dtemp_13-dtemp_14);
  }

  return jacobi;
}

double eval_grad_jacobi(double in_r, int in_alpha, int in_beta, int in_mode)
{
  double grad_jacobi;

  if (in_mode==0){
    grad_jacobi = 0.;
  } else {
    grad_jacobi=sqrt(1.0*in_mode*(in_mode+in_alpha+in_beta+1))*eval_jacobi(in_r,in_alpha+1,in_beta+1,in_mode-1);
  }

  return grad_jacobi;
}

double eval_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order)
{
  double dubiner_basis_2d;

  int n_dof = ((in_basis_order+1)*(in_basis_order+2))/2;

  if(in_mode<n_dof) {
    int i,j,k;
    int mode;
    double jacobi_0, jacobi_1;
    point ab;

    ab = rs_to_ab(in_rs);

    mode = 0;
    for (k=0;k<in_basis_order+1;k++) {
      for (j=0;j<k+1;j++) {
        i = k-j;

        if(mode==in_mode) {
          // found the correct mode
          jacobi_0=eval_jacobi(ab.x,0,0,i);
          jacobi_1=eval_jacobi(ab.y,(2*i)+1,0,j);
          dubiner_basis_2d=sqrt(2.0)*jacobi_0*jacobi_1*pow(1.0-ab.y,i);
        }

        mode++;
      }
    }
  } else {
    FatalError("Invalid mode when evaluating Dubiner basis.");
  }

  return dubiner_basis_2d;
}

double eval_dr_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order)
{
  double dr_dubiner_basis_2d;

  int n_dof=((in_basis_order+1)*(in_basis_order+2))/2;

  if(in_mode<n_dof) {
    int i,j,k;
    int mode;
    double jacobi_0, jacobi_1;
    point ab;

    ab = rs_to_ab(in_rs);

    mode = 0;
    for (k=0;k<in_basis_order+1;k++) {
      for (j=0;j<k+1;j++) {
        i = k-j;

        if(mode==in_mode) {
          // found the correct mode
          jacobi_0=eval_grad_jacobi(ab.x,0,0,i);
          jacobi_1=eval_jacobi(ab.y,(2*i)+1,0,j);

          if(i==0) {
            // to avoid singularity
            //dr_dubiner_basis_2d=sqrt(2.0)*jacobi_0*jacobi_1;
            dr_dubiner_basis_2d=0.;
          } else {
            dr_dubiner_basis_2d=2.0*sqrt(2.0)*jacobi_0*jacobi_1*pow(1.0-ab.y,i-1);
          }
        }

        mode++;
      }
    }
  } else {
    FatalError("Invalid mode when evaluating basis.");
  }

  return dr_dubiner_basis_2d;
}

// helper method to evaluate d/ds of scalar dubiner basis

double eval_ds_dubiner_basis_2d(point &in_rs, int in_mode, int in_basis_order)
{
  double ds_dubiner_basis_2d;

  int n_dof=((in_basis_order+1)*(in_basis_order+2))/2;

  if(in_mode<n_dof)
  {
    int i,j,k;
    int mode;
    double jacobi_0, jacobi_1, jacobi_2, jacobi_3, jacobi_4;
    point ab;

    ab = rs_to_ab(in_rs);


    mode = 0;
    for (k=0;k<in_basis_order+1;k++) {
      for (j=0;j<k+1;j++) {
        i = k-j;

        if(mode==in_mode)  {
          // found the correct mode
          jacobi_0 = eval_grad_jacobi(ab.x,0,0,i);
          jacobi_1 = eval_jacobi(ab.y,(2*i)+1,0,j);

          jacobi_2 = eval_jacobi(ab.x,0,0,i);
          jacobi_3 = eval_grad_jacobi(ab.y,(2*i)+1,0,j)*pow(1.0-ab.y,i);
          jacobi_4 = eval_jacobi(ab.y,(2*i)+1,0,j)*i*pow(1.0-ab.y,i-1);

          if (i==0) {
            // to avoid singularity
            ds_dubiner_basis_2d = sqrt(2.0)*(jacobi_2*jacobi_3);
          } else {
            ds_dubiner_basis_2d = sqrt(2.0)*((jacobi_0*jacobi_1*pow(1.0-ab.y,i-1)*(1.0+ab.x))+(jacobi_2*(jacobi_3-jacobi_4)));
          }
        }

        mode++;
      }
    }
  }
  else
  {
    FatalError("ERROR: Invalid mode when evaluating basis.");
  }

  return ds_dubiner_basis_2d;
}

double eval_gamma(int in_n)
{
  int i;
  double gamma_val;

  if(in_n==1)
    gamma_val=1;
  else
  {
    gamma_val=in_n-1;

    for(i=0; i<in_n-2; i++) {
      gamma_val=gamma_val*(in_n-2-i);
    }
  }

  return gamma_val;
}

point rs_to_ab(point &rs)
{
  point ab;

  if (rs.x == 1.0) { // to avoid singularity
    ab.x = -1.0;
  } else {
    ab.y = (2.0*((1.0+rs.x)/(1.0-rs.y))) - 1.0;
  }

  ab.y = rs.y;

  return ab;
}
