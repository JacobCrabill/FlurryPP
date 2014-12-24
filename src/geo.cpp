/*!
 * \file geo.cpp
 * \brief Class for handling geometry setup & modification
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

#include "../include/geo.hpp"


geo::setup(input* _params)
{
  params = _params;

  switch(params->mesh_type) {
  case (READ_MESH):
    readGmsh(params->mesh_file_name);

  case (CREATE_MESH):
    createMesh();

  default:
    FatalError("Mesh type not recognized.")
  }
}


geo::processConnectivity()
{
  /* --- Setup Edges --- */
  vector<int*> e2v;
  int iep1, ne = 0;
  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2nv(e); ie++) {
      int* edge = new int(2);
      iep1 = (ie+1)%c2nv(e);
      edge[0] = c2v[e][ie];
      edge[1] = c2v[e][iep1];
      e2v.insert(edge)
    }
  }

}

geo::setupElesFaces(solver *Solver)
{
  if (nEles<=0) FatalError("Cannot setup elements array - nEles = 0");

  // Temporary vectors that will get assigned to Solver
  vector<ele> eles(nEles);
  vector<face> faces(nFaces);

  // Setup the elements
  ele tmp;
  for (int e=0; e<nEles; e++) {
    tmp.ID = e;

    // Shape [mesh] nodes
    tmp.nodeID.resize(c2nv[e]);
    tmp.nodes.resize(c2nv[e]);
    for (int iv=0; iv<c2nv[e]; iv++) {
      tmp.nodeID[iv] = c2v[e][iv];
      tmp.nodes[iv] = xv[c2v[e][iv]];
    }

    // Global face IDs
    tmp.faceID.resize(c2ne[e]);
    for (int ie=0; ie<c2ne[e]; ie++) {
      tmp.faceID[ie] = c2e[e][ie];
    }

    eles[e] = tmp;
  }

  // Final Step - Assign eles, faces to Solver
  Solver->eles = eles;
  Solver->faces = faces;
}

geo::createMesh()
{
  int i, j, nx, ny;
  int ne, nv;
  double xmin, xmax, ymin, ymax, dx, dy;
  vector<int> c2v_tmp;

  nx = params->nx;
  ny = params->ny;

  xmin = params->xmin;
  xmax = params->xmax;
  ymin = params->ymin;
  ymax = params->ymax;

  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;

  nEles = nx*ny;
  nVerts = (nx+1)*(ny+1);

  xv.resize(nVerts);
  c2v_tmp.resize(4);

  /* --- Setup Vertices --- */
  nv = 0;
  point pt;
  for (i=0; i<nx+1; i++) {
    for (j=0; j<ny+1; j++) {
      pt.x = xmin + i*dx;
      pt.y = ymin + j*dy;
      xv[nv] = pt;
      nv++;
    }
  }

  /* --- Setup Elements --- */
  ne = 0;
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      c2v_tmp[0] = j*(nx+1) + i;
      c2v_tmp[1] = j*(nx+1) + i + 1;
      c2v_tmp[2] = (j+1)*(nx+1) + i + 1;
      c2v_tmp[3] = (j+1)*(nx+1) + i;
      c2v.insertRow(c2v_tmp);
      ne++;
    }
  }
}
