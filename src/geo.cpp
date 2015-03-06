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
      if (c2v[e][ie] < c2v[e][iep1]) {
      edge[0] = c2v[e][ie];
      edge[1] = c2v[e][iep1];
      } else {
        edge[0] = c2v[e][iep1];
        edge[1] = c2v[e][ie];
      }
      e2v1.insertRow(edge);
    }
  }

  /* --- Just for nDims==2 ---
     Get just the unique edges
     iE is of length [original e2v] with range [final e2v] */
  vector<int> iE;
  e2v1.unique(e2v,iE);
  nEdges = e2v.getDim0();

  /* --- Determine interior vs. boundary edges/faces --- */

  // Flag for whether global face ID corresponds to interior or boundary face
  isBnd.assign(nEdges,0);

  nFaces = 0;
  nBndFaces = 0;

  vector<int> ie;
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
        nFaces++;
      }
      else if (ie.size()==1) {
        // Boundary Edge
        bndEdges.push_back(iE[i]); //[nBndEdges] = iE[i];
        isBnd[iE[i]] = 1;
        nBndFaces++;
      }

      // Mark edges as completed
      vecAssign(iE,ie,-1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bcType.resize(nBndFaces);


  /* --- Setup Edge Normals? --- */

  /* --- Setup Cell-To-Edge, Edge-To-Cell --- */
  c2e.setup(nEles,getMax(c2ne));
  c2b.setup(nEles,getMax(c2ne));
  c2b.initializeToZero();
  e2c.setup(nEdges,2);
  e2c.initializeToValue(-1);

  vector<int> edge(2), ie1(2), col2;
  int ie0, ie2;
  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2ne[ic]; j++) {
      int jp1 = (j+1)%(c2ne[ic]);

      // Store edges consistently to allow matching of duplicates
      if (c2v[ic][j] < c2v[ic][jp1]) {
      edge[0] = c2v[ic][j];
      edge[1] = c2v[ic][jp1];
      } else {
        edge[0] = c2v[ic][jp1];
        edge[1] = c2v[ic][j];
      }

      ie1 = findEq(e2v.getCol(0),edge[0]);
      col2 = (e2v.getRows(ie1)).getCol(1);
      ie2 = findFirst(col2,edge[1]);
      ie0 = ie1[ie2];

      // Find ID of face within type-specific array
      if (isBnd[ie0]) {
        c2e[ic][j] = ie0; //findFirst(bndEdges,ie0);
        c2b[ic][j] = 1;
      }else{
        c2e[ic][j] = ie0; //findFirst(intEdges,ie0);
        c2b[ic][j] = 0;
      }

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

void geo::setupElesFaces(vector<ele> &eles, vector<face> &faces, vector<bound> bounds)
{
  if (nEles<=0) FatalError("Cannot setup elements array - nEles = 0");

  eles.resize(nEles);
  faces.resize(nFaces);
  bounds.resize(nBndFaces);

  // Setup the elements
  int ic = 0;
  for (auto& e:eles) {
    e.ID = ic;
    e.eType = ctype[ic];
    e.nNodes = c2nv[ic];

    // Shape [mesh] nodes
    e.nodeID.resize(c2nv[ic]);
    e.nodes.resize(c2nv[ic]);
    for (int iv=0; iv<c2nv[ic]; iv++) {
      e.nodeID[iv] = c2v[ic][iv];
      e.nodes[iv] = xv[c2v[ic][iv]];
    }

    // Global face IDs for internal & boundary faces
    e.faceID.resize(c2ne[ic]);
    e.bndFace.resize(c2ne[ic]);
    for (int k=0; k<c2ne[ic]; k++) {
      e.bndFace[k] = c2b[ic][k];
      e.faceID[k] = c2e[ic][k];
    }

    e.setup(params,this);

    ic++;
  }

  /* --- Setup the faces --- */

  vector<int> tmpEdges;
  int fid1, fid2;
  int ie, i = 0;

  // Internal Faces
  for (auto& F:faces) {
    // Find global face ID of current interior face
    ie = intEdges[i];
    ic = e2c[ie][0];
    // Find local face ID of global face within first element [on left]
    tmpEdges = c2e[ic];
    fid1 = findFirst(tmpEdges,ie);
    F.params = params;
    if (e2c[ie][1] == -1) {
      FatalError("Interior edge does not have a right element assigned.");
    }else{
      tmpEdges = c2e[e2c[ie][1]];
      fid2 = findFirst(tmpEdges,ie);
      F.setupFace(&eles[e2c[ie][0]],&eles[e2c[ie][1]],fid1,fid2,ie);
    }

    i++;
  }

  // Boundary Faces
  i = 0;
  for (auto& B:bounds) {
    // Find global face ID of current boundary face
    ie = bndEdges[i];
    ic = e2c[ie][0];
    // Find local face ID of global face within element
    tmpEdges = c2e[ic];
    fid1 = findFirst(tmpEdges,ie);
    B.params = params;
    if (e2c[ie][1] != -1) {
      FatalError("Boundary edge has a right element assigned.");
    }else{
      B.setupBound(&eles[e2c[ie][0]],fid1,bcType[i],ie);
    }

    i++;
  }
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
