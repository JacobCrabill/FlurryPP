/*!
 * \file superMesh.cpp
 * \brief superMesh class definition
 *
 * Creates a local supermesh for an element in one grid from
 * elements in another grid (See Farrell and Maddison, 2010)
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill
 *
 */

#include "superMesh.hpp"

#include "funcs.hpp"

superMesh::superMesh()
{

}

superMesh::~superMesh()
{

}

void superMesh::setup(geo* _gridT, geo* _gridD, int _targetCell, int _order)
{
  gridT = _gridT;
  gridD = _gridD;

  targetCell = _targetCell;
  order = _order;

  buildSuperMesh();
}

void superMesh::buildSuperMesh(void)
{

}

Array<double,3> superMesh::splitHexIntoTet(matrix<double> &hexNodes)
{
  Array<double,3> tetNodes(5,4,3);

  vector<vector<int>> ind = {{0,1,4,3},{2,1,6,3},{5,1,6,4},{7,3,4,6},{1,3,6,4}};

  for (int i=0; i<5; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<3; k++)
        tetNodes(i,j,k) = hexNodes(ind[i][j],k);

  return tetNodes;
}

Array<double,3> superMesh::clipTet(matrix<double> &tetNodes, matrix<double> &clipFace, Vec3 &norm)
{

}

double superMesh::integrate(vector<double> &data)
{
  if (data.size() != nQpts) FatalError("To integrate over supermesh, data must lie at its quadrature nodes.");

  double val = 0;
  for (int i=0; i<nTets; i++) {
    for (int j=0; j<nQpts_tet; j++) {

    }
  }
}
