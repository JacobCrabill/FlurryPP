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


double Lagrange(vector<double> &x_lag, double &y, int &mode)
{
  int i;
  double lag = 1.0;

  for(i=0; i<x_lag.size(); i++) {
    if(i!=mode) {
      lag = lag*((yr-x_lag[i])/(x_lag(mode)-x_lag[i]));
    }
  }

  return lag;
}


double dLagrange(vector<double> &x_lag, double &y, int &mode)
{
  int i, j;
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
          dLag_den = dLag_den*(x_lag(mode)-x_lag[j]);
        }
      }

      dLag = dLag+(dLag_num/dLag_den);
    }
  }

  return dLag;
}

double ddLagrange(vector<double> &x_lag, double &y, int &mode)
{
  int i, j, k;
  double ddLag, ddLag_num, ddLag_den;

  ddLag = 0.0;

  for (i=0; i<x_lag.size(); i++) {
    if (i!=mode) {

      for (j=0; j<x_lag.size(); j++) {
        if (j!=mode) {
          ddLag_num = 1.0;
          ddLag_den = 1.0;

          for (k=0; k<x_lag.size(); k++) {
            if(k!=in_mode) {
              if(k!=i && k!=j) {
                dtemp_1 = dtemp_1*(y-x_lag[k]);
              }

              dtemp_2 = dtemp_2*(x_lag[mode]-x_lag[k]);
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
  out_shape.resize(4);

  out_shape[0] = 0.25*(1-in_rs.x)*(1-in_rs.y);
  out_shape[1] = 0.25*(1+in_rs.x)*(1-in_rs.y);
  out_shape[2] = 0.25*(1+in_rs.x)*(1+in_rs.y);
  out_shape[3] = 0.25*(1-in_rs.x)*(1+in_rs.y);
}

void dshape_quad(point &in_rs, vector<vector<double>> &out_dshape)
{
  out_dshape.resize(4);
  for (int i=0; i<4; i++)
    out_dshape[i].resize(2);

  out_dshape[0][0] = -0.25*(1-in_rs.y);
  out_dshape[1][0] =  0.25*(1-in_rs.y);
  out_dshape[2][0] =  0.25*(1+in_rs.y;
  out_dshape[3][0] = -0.25*(1+in_rs.y);

  out_dshape[0][1] = -0.25*(1-in_rs.x);
  out_dshape[1][1] = -0.25*(1+in_rs.x);
  out_dshape[2][1] =  0.25*(1+in_rs.x);
  out_dshape[3][1] =  0.25*(1-in_rs.x);
}
