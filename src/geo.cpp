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

geo::geo()
{

}

void geo::setup(input* params)
{
  this->params = params;

  switch(params->mesh_type) {
    case (READ_MESH):
      readGmsh(params->meshFileName);
      break;

    case (CREATE_MESH):
      createMesh();
      break;

    default:
      FatalError("Mesh type not recognized.");
  }

  processConnectivity();
}


void geo::processConnectivity()
{
  /* --- Setup Edges --- */
  matrix<int> e2v1, e2v;
  int iep1;
  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2nv[e]; ie++) {
      vector<int> edge(2);
      iep1 = (ie+1)%c2nv[e];
      edge[0] = c2v[e][ie];
      edge[1] = c2v[e][iep1];
      e2v1.insertRow(edge);
    }
  }

  /* --- Just for nDims==2 ---
     Get just the unique edges
     iE is of length [original e2v] with range [final e2v] */
  vector<int> iE;
  e2v1.unique(e2v,iE);
  nEdges = e2v.getDim0();

  if (nDims == 2) nFaces = nEdges;

  vector<int> ie;
  vector<int> intEdges, bndEdges; // Setup something more permanent...
  for (uint i=0; i<iE.size(); i++) {
    if (iE[i]!=-1) {
      ie = findEq(iE,iE[i]);
      if (ie.size()>2) {
        string errMsg = "More than 2 cells for edge " + to_string(i);
        FatalError(errMsg.c_str());
      }
      else if (ie.size()==2) {
        // Internal Edge which has not yet been added
        intEdges.push_back(iE[i]); //[nIntEdges] = iE[i];
      }
      else if (ie.size()==1) {
        // Boundary Edge
        bndEdges.push_back(iE[i]); //[nBndEdges] = iE[i];
      }

      // Mark edges as completed
      vecAssign(iE,ie,-1);
    }
  }

  /* --- Setup Edge Normals? --- */

  /* --- Setup Cell-To-Edge, Edge-To-Cell --- */
  c2e.setup(nEles,getMax(c2ne));
  e2c.setup(nEdges,2);
  e2c.initializeToValue(-1);

  vector<int> edge(2), ie1(2), col2;
  int ie0, ie2;
  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2ne[ic]; j++) {
      int jp1 = (j+1)%(c2ne[ic]);
      edge[0] = c2v[ic][j];
      edge[1] = c2v[ic][jp1];

      ie1 = findEq(e2v.getCol(0),edge[0]);
      col2 = (e2v.getRows(ie1)).getCol(1);
      ie2 = findFirst(col2,edge[1]);
      ie0 = ie1[ie2];

      c2e[ic][j] = ie0;
      if (e2c[ie0][0] == -1) {
        // No cell yet assigned to edge; put on left
        e2c[ie0][0] = ic;
      }else{
        // Put cell on right
        e2c[ie0][1] = ic;
      }
    }
  }
}

void geo::setupElesFaces(solver *Solver)
{
  if (nEles<=0) FatalError("Cannot setup elements array - nEles = 0");

  // Temporary vectors that will get assigned to Solver
  vector<ele> eles(nEles);
  vector<face> faces(nFaces);

  // Setup the elements
  int ie = 0;
  for (auto& e:eles) {
    e.ID = ie;
    e.eType = ctype[ie];
    e.nNodes = c2nv[ie];

    // Shape [mesh] nodes
    e.nodeID.resize(c2nv[ie]);
    e.nodes.resize(c2nv[ie]);
    for (int iv=0; iv<c2nv[ie]; iv++) {
      e.nodeID[iv] = c2v[ie][iv];
      e.nodes[iv] = xv[c2v[ie][iv]];
    }

    // Global face IDs
    e.faceID.resize(c2ne[ie]);
    for (int k=0; k<c2ne[ie]; k++) {
      e.faceID[k] = c2e[ie][k];
    }

    e.setup(params,this);

    ie++;
  }

  // Setup the faces
  vector<int> tmpEdges;
  int fid1, fid2;
  int i = 0;
  for (auto& F:faces) {
    // Find local face ID of global face
    tmpEdges = c2e[e2c[i][0]];
    fid1 = findFirst(tmpEdges,i);
    F.params = params;
    if (e2c[i][1] == -1) {
      // Boundary Edge
      // ** still need to create a boundary face class (and an interior face class?) **
    }else{
      // Interior Edge
      tmpEdges = c2e[e2c[i][1]];
      fid2 = findFirst(tmpEdges,i);
      F.setupFace(&eles[e2c[i][0]],&eles[e2c[i][1]],fid1,fid2,i);
    }

    i++;
  }

  // Final Step - Assign eles, faces to Solver
  Solver->eles = eles;
  Solver->faces = faces;
}

void geo::readGmsh(string fileName)
{

}

void geo::createMesh()
{
  int i, j, nx, ny;
  int ne, nv;
  double xmin, xmax, ymin, ymax, dx, dy;
  vector<int> c2v_tmp;

  nx = params->nx;
  ny = params->ny;
  nDims = params->nDims; // Since I may implement createMesh for 3D later

  xmin = params->xmin;
  xmax = params->xmax;
  ymin = params->ymin;
  ymax = params->ymax;

  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;

  nEles = nx*ny;
  nVerts = (nx+1)*(ny+1);

  c2nv.assign(nEles,4);
  c2ne.assign(nEles,4);
  ctype.assign(nEles,QUAD); // Add hex later

  xv.resize(nVerts);
  c2v_tmp.resize(4);

  /* --- Setup Vertices --- */
  nv = 0;
  point pt;
  for (i=0; i<ny+1; i++) {
    for (j=0; j<nx+1; j++) {
      pt.x = xmin + j*dx;
      pt.y = ymin + i*dy;
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

#include "../include/geo.inl"
