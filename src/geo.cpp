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

  processPeriodicBoundaries();
}


void geo::processConnectivity()
{
  /* --- Setup Edges --- */
  matrix<int> e2v1, e2v;
  vector<int> edge(2);
  int iep1;

  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2nv[e]; ie++) {
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

  /* --- Determine interior vs. boundary edges/faces --- */

  // Flag for whether global face ID corresponds to interior or boundary face
  isBnd.assign(nEdges,0);

  nFaces = 0;
  nBndEdges = 0;

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
        intEdges.push_back(iE[i]);
        nFaces++;
      }
      else if (ie.size()==1) {
        // Boundary Edge
        bndEdges.push_back(iE[i]);
        isBnd[iE[i]] = true;
        nBndEdges++;
      }

      // Mark edges as completed
      vecAssign(iE,ie,-1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bcType.assign(nBndEdges,-1);
  for (int i=0; i<nBndEdges; i++) {
    int iv1 = e2v[bndEdges[i]][0];
    int iv2 = e2v[bndEdges[i]][1];
    for (int bnd=0; bnd<nBounds; bnd++) {
      if (findFirst(bndPts[bnd],iv1)!=-1 && findFirst(bndPts[bnd],iv2)!=-1) {
        // The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bndFaces[bnd].insertRow(e2v[bndEdges[i]]);
        break;
      }
    }
  }

  /* --- Setup Cell-To-Edge, Edge-To-Cell --- */

  c2e.setup(nEles,getMax(c2ne));
  c2b.setup(nEles,getMax(c2ne));
  c2b.initializeToZero();
  e2c.setup(nEdges,2);
  e2c.initializeToValue(-1);

  vector<int> ie1(2), col2;
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

      // Find ID of face within type-specific array
      if (isBnd[ie0]) {
        c2e[ic][j] = findFirst(bndEdges,ie0);
        c2b[ic][j] = true;
      }else{
        c2e[ic][j] = findFirst(intEdges,ie0);
        c2b[ic][j] = false;
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
  bounds.resize(nBndEdges);

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
    ie = e2c[intEdges[i]][0];
    // Find local face ID of global face within first element [on left]
    tmpEdges = c2e[ie];
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
    ie = e2c[bndEdges[i]][0];
    // Find local face ID of global face within element
    tmpEdges = c2e[ie];
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

  params->periodicDX = xmax-xmin;
  params->periodicDY = ymax-ymin;

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
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      c2v_tmp[0] = j*(nx+1) + i;
      c2v_tmp[1] = j*(nx+1) + i + 1;
      c2v_tmp[2] = (j+1)*(nx+1) + i + 1;
      c2v_tmp[3] = (j+1)*(nx+1) + i;
      c2v.insertRow(c2v_tmp);
    }
  }

  /* --- Setup Boundaries --- */
  int bc = PERIODIC;
  nBounds = 1;
  bcList.push_back(bc);
  bndFaces.resize(1);
  bndPts.setup(1,(nx+1)*(ny+1));
  bndFaces[0].setup(2*nx+2*ny,2); // nBndFaces x 2-pts-per-edge
  ne = 0;

  // Bottom Edge Faces
  for (int ix=0; ix<nx; ix++) {
    bndFaces[0][ne][0] = ix;
    bndFaces[0][ne][1] = ix+1;
    bndPts[0][2*ne]   = bndFaces[0][ne][0];
    bndPts[0][2*ne+1] = bndFaces[0][ne][1];
    ne++;
  }

  // Top Edge Faces
  for (int ix=0; ix<nx; ix++) {
    bndFaces[0][ne][1] = (nx+1)*ny + ix;
    bndFaces[0][ne][0] = (nx+1)*ny + ix+1;
    bndPts[0][2*ne]   = bndFaces[0][ne][0];
    bndPts[0][2*ne+1] = bndFaces[0][ne][1];
    ne++;
  }

  // Left Edge Faces
  for (int iy=0; iy<ny; iy++) {
    bndFaces[0][ne][1] = iy*(nx+1);
    bndFaces[0][ne][0] = (iy+1)*(nx+1);
    bndPts[0][2*ne]   = bndFaces[0][ne][0];
    bndPts[0][2*ne+1] = bndFaces[0][ne][1];
    ne++;
  }

  // Right Edge Faces
  for (int iy=0; iy<ny; iy++) {
    bndFaces[0][ne][0] = iy*(nx+1) + nx;
    bndFaces[0][ne][1] = (iy+1)*(nx+1) + nx;
    bndPts[0][2*ne]   = bndFaces[0][ne][0];
    bndPts[0][2*ne+1] = bndFaces[0][ne][1];
    ne++;
  }

  // Remove duplicates in bndPts
  std::sort(bndPts[0].begin(), bndPts[0].end());
  vector<int>::iterator it;
  it = std::unique (bndPts[0].begin(), bndPts[0].end());
  bndPts[0].resize(std::distance(bndPts[0].begin(), it));

  // Default created mesh will be a rectangle with periodic boundaries
  // Top, Right edges will be replace by equivalent Bottom, Left edges
//  int ic;
//  for (int iy=0; iy<ny-1; iy++) {
//    for (int ix=0; ix<nx-1; ix++) {
//      ic = ix + ny*iy;
//      for (int ie=0; ie<4; ie++) {
//        iep1 = (ie+1)%4;
//        edge[0] = c2v[ic][ie];
//        edge[1] = c2v[ic][iep1];
//        e2v1.insertRow(edge);
//      }
//    }
//    // Match right edges of mesh @ row iy to left edge
//    ic = nx + ny*iy;
//    edge[0] = c2v[ic][0];        // Bottom edge of current cell
//    edge[1] = c2v[ic][1];
//    e2v1.insertRow(edge);
//    edge[0] = c2v[ic-(nx-1)][3]; // Left edge of left-most cell in row
//    edge[1] = c2v[ic-(nx-1)][0];
//    e2v1.insertRow(edge);
//    edge[0] = c2v[ic][2];        // Top edge of current cell
//    edge[1] = c2v[ic][3];
//    e2v1.insertRow(edge);
//  }

//  // Match top edges of mesh to bottom edges
//  ic = ny*(ny-1);
//  edge[0] = c2v[ic][3];
//  edge[1] = c2v[ic][0];
//  e2v1.insertRow(edge);
//  for (int ix=0; ix<nx-1; ix++) {
//    ic = ix + ny*(ny-1);
//    edge[0] = c2v[ic][1]; // Right edge
//    edge[1] = c2v[ic][2];
//    e2v1.insertRow(edge);
//    edge[0] = c2v[ix][0]; // Bottom edge of cell at bottom of mesh
//    edge[1] = c2v[ix][1];
//    e2v1.insertRow(edge);
//  }

//  ic = nx*(ny-1); // Top-left cell
//  edge[0] = c2v[ic][3];   // Final cell, right egde
//  edge[1] = c2v[ic][0];
//  e2v1.insertRow(edge);
//  edge[0] = nx-1; // Final cell, top edge
//  edge[1] = nx;
//  e2v1.insertRow(edge);
}

void geo::processPeriodicBoundaries(void)
{
  uint nPeriodic;
  vector<int> iPeriodic(0);

  for (int i=0; i<nBndEdges; i++) {
    if (bcType[i] == PERIODIC) {
      iPeriodic.push_back(i);
    }
  }

  nPeriodic = iPeriodic.size();
  if (nPeriodic == 0) return;
  if (nPeriodic%2 != 0) FatalError("Expecting even number of periodic faces; have odd number.");

  for (auto& i:iPeriodic) {
    if (bndEdges[i]==-1) continue;
    for (auto& j:iPeriodic) {
      if (i==j || bndEdges[j]==-1) continue;
      if (checkPeriodicFaces(e2v[bndEdges[i]],e2v[bndEdges[j]])) {
        /* --- Match found - now take care of transfer from boundary -> internal --- */

        if (i>j) FatalError("How did this happen?!");

        // Transfer combined edge from boundary to internal list
        intEdges.push_back(bndEdges[i]);

        // Flag global edge IDs as internal edges
        isBnd[bndEdges[i]] = false;
        isBnd[bndEdges[j]] = false;

        // Fix e2c - add right cell to combined edge, make left cell = -1 in 'deleted' edge
        e2c[bndEdges[i]][1] = e2c[bndEdges[j]][0];
        e2c[bndEdges[j]][0] = -1;

        // Fix c2e - replace 'deleted' edge from right cell with combined edge
        int fID = findFirst(c2e[e2c[i][1]],bndEdges[j]);
        c2e[e2c[i][1]][fID] = bndEdges[i];

        // Fix c2b - set element-local face to be internal face
        c2b[e2c[i][1]][fID] = false;

        // Flag edges as gone in boundary edges list
        bndEdges[i] = -1;
        bndEdges[j] = -1;
      }
    }
  }
}

bool geo::checkPeriodicFaces(vector<int> edge1, vector<int> edge2)
{
  double x11, x12, y11, y12, x21, x22, y21, y22;
  x11 = xv[edge1[0]][0];  y11 = xv[edge1[0]][1];
  x12 = xv[edge1[1]][0];  y12 = xv[edge1[1]][1];
  x21 = xv[edge2[0]][0];  y21 = xv[edge2[0]][1];
  x22 = xv[edge2[1]][0];  y22 = xv[edge2[1]][1];

  double tol = params->periodicTol;
  double dx = params->periodicDX;
  double dy = params->periodicDY;

  if ( abs(abs(x21-x11)-dx)<tol && abs(abs(x22-x12)-dx)<tol && abs(y21-y11)<tol && abs(y22-y12)<tol ) {
    // Faces match up across x-direction, with [0]->[0] and [1]->[1] in each edge
    return true;
  }
  else if ( abs(abs(x22-x11)-dx)<tol && abs(abs(x21-x12)-dx)<tol && abs(y22-y11)<tol && abs(y21-y12)<tol ) {
    // Faces match up across x-direction, with [0]->[1] and [1]->[0] in each edge
    return true;
  }
  else if ( abs(abs(y21-y11)-dy)<tol && abs(abs(y22-y12)-dy)<tol && abs(x21-x11)<tol && abs(x22-x12)<tol ) {
    // Faces match up across y-direction, with [0]->[0] and [1]->[1] in each edge
    return true;
  }
  else if ( abs(abs(y22-y11)-dy)<tol && abs(abs(y21-y12)-dy)<tol && abs(x22-x11)<tol && abs(x21-x12)<tol ) {
    // Faces match up across y-direction, with [0]->[1] and [1]->[0] in each edge
    return true;
  }
  else {
    return false;
  }
}

#include "../include/geo.inl"
