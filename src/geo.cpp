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

#include <algorithm>
#include <cstdlib>

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
  matrix<int> e2v1;
  vector<int> edge(2);

  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2nv[e]; ie++) {
      int iep1 = (ie+1)%c2nv[e];
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

  /* --- Just for nDims==2: Get just the unique edges --- */
  // iE is of length [original e2v1] with range [final e2v]
  // The number of times an edge appears in iE is equal to
  // the number of cells that edge touches
  vector<int> iE;
  e2v1.unique(e2v,iE);
  nEdges = e2v.getDim0();

  /* --- Determine interior vs. boundary edges/faces --- */

  // Flag for whether global face ID corresponds to interior or boundary face
  isBnd.assign(nEdges,0);

  nFaces = 0;
  nBndEdges = 0;

  for (uint i=0; i<iE.size(); i++) {
    if (iE[i]!=-1) {
      vector<int> ie = findEq(iE,iE[i]);
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
  // Since bndFaces was setup during createMesh, and will be created
  // here anyways, clear bndFaces & re-setup
  bndFaces.clear();
  bndFaces.resize(nBounds);
  bcType.assign(nBndEdges,-1);
  for (int i=0; i<nBndEdges; i++) {
    int iv1 = e2v[bndEdges[i]][0];
    int iv2 = e2v[bndEdges[i]][1];
    for (int bnd=0; bnd<nBounds; bnd++) {
      if (findFirst(bndPts[bnd],iv1,bndPts.dim1)!=-1 && findFirst(bndPts[bnd],iv2,bndPts.dim1)!=-1) {
        // The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bndFaces[bnd].insertRow(e2v[bndEdges[i]],-1,e2v.dim1);
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

      vector<int> ie1 = findEq(e2v.getCol(0),edge[0]);
      vector<int> col2 = (e2v.getRows(ie1)).getCol(1);
      int ie2 = findFirst(col2,edge[1]);
      int ie0 = ie1[ie2];

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

void geo::setupElesFaces(vector<ele> &eles, vector<face> &faces, vector<bound> &bounds)
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
  int i = 0;

  // Internal Faces
  for (auto& F:faces) {
    // Find global face ID of current interior face
    int ie = intEdges[i];
    ic = e2c[ie][0];
    // Find local face ID of global face within first element [on left]
    tmpEdges.assign(c2e[ic],c2e[ic]+c2ne[ic]);
    int fid1 = findFirst(tmpEdges,ie);
    F.params = params;
    if (e2c[ie][1] == -1) {
      FatalError("Interior edge does not have a right element assigned.");
    }else{
      ic = e2c[ie][1];
      tmpEdges.assign(c2e[ic], c2e[ic]+c2ne[ic]);
      //tmpEdges = c2e[e2c[ie][1]]; // previous row-as-vector form
      int fid2 = findFirst(tmpEdges,ie);
      F.setupFace(&eles[e2c[ie][0]],&eles[e2c[ie][1]],fid1,fid2,ie);
    }

    i++;
  }

  // Boundary Faces
  i = 0;
  for (auto& B:bounds) {
    // Find global face ID of current boundary face
    int ie = bndEdges[i];
    ic = e2c[ie][0];
    // Find local face ID of global face within element
    tmpEdges.assign(c2e[ic],c2e[ic]+c2ne[ic]);
    int fid1 = findFirst(tmpEdges,ie);
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
  int nx = params->nx;
  int ny = params->ny;
  nDims = params->nDims; // Since I may implement createMesh for 3D later

  double xmin = params->xmin;
  double xmax = params->xmax;
  double ymin = params->ymin;
  double ymax = params->ymax;

  double dx = (xmax-xmin)/nx;
  double dy = (ymax-ymin)/ny;

  params->periodicDX = xmax-xmin;
  params->periodicDY = ymax-ymin;

  nEles = nx*ny;
  nVerts = (nx+1)*(ny+1);

  c2nv.assign(nEles,4);
  c2ne.assign(nEles,4);
  ctype.assign(nEles,QUAD); // Add hex later

  xv.resize(nVerts);
  vector<int> c2v_tmp(4,0);

  /* --- Setup Vertices --- */
  int nv = 0;
  point pt;
  for (int i=0; i<ny+1; i++) {
    for (int j=0; j<nx+1; j++) {
      pt.x = xmin + j*dx;
      pt.y = ymin + i*dy;
      xv[nv] = pt;
      nv++;
    }
  }

//  // Add a random perturbation of +/- dx/4 to interior points
//  for (i=1; i<ny; i++) {
//    for (j=1; j<nx; j++) {
//      pt.x = (dx/4)*(-1+((double)(rand()%1000))/500.);
//      pt.y = (dy/4)*(-1+((double)(rand()%1000))/500.);
//      xv[i*(ny+1)+j].x += pt.x;
//      xv[i*(ny+1)+j].y += pt.y;
//    }
//  }

  /* --- Setup Elements --- */
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      c2v_tmp[0] = j*(nx+1) + i;
      c2v_tmp[1] = j*(nx+1) + i + 1;
      c2v_tmp[2] = (j+1)*(nx+1) + i + 1;
      c2v_tmp[3] = (j+1)*(nx+1) + i;
      c2v.insertRow(c2v_tmp);
    }
  }

  /* --- Setup Boundaries --- */
  // List of all boundary conditions being used (bcNum maps string->int)
  bcList.push_back(bcNum[params->create_bcBottom]);
  bcList.push_back(bcNum[params->create_bcRight]);
  bcList.push_back(bcNum[params->create_bcTop]);
  bcList.push_back(bcNum[params->create_bcLeft]);

  // Sort the list & remove any duplicates
  std::sort(bcList.begin(), bcList.end());
  vector<int>::iterator vIt = std::unique(bcList.begin(), bcList.end());
  nBounds = std::distance(bcList.begin(), vIt);     // will I need both an nBounds (i.e., in mesh) and an nBC's (current nBounds)?
  bcList.resize(nBounds);

  // Setup a map so we know where each BC# is inside of bcList
  map<int,int> bc2bcList;

  // Setup boundary connectivity storage
  nBndFaces.assign(nBounds,0);
  bndFaces.resize(nBounds);
  bndPts.setup(nBounds,4*nx+4*ny); //(nx+1)*(ny+1));
  nBndPts.resize(nBounds);
  for (int i=0; i<nBounds; i++) {
    bc2bcList[bcList[i]] = i;
    bndFaces[i].setup(2*nx+2*ny,2); // max nBndFaces x 2-pts-per-edge
  }

  // Bottom Edge Faces
  int ib = bc2bcList[bcNum[params->create_bcBottom]];
  int ne = nBndFaces[ib];
  for (int ix=0; ix<nx; ix++) {
    bndFaces[ib][ne][0] = ix;
    bndFaces[ib][ne][1] = ix+1;
    bndPts[ib][2*ne]   = bndFaces[ib][ne][0];
    bndPts[ib][2*ne+1] = bndFaces[ib][ne][1];
    ne++;
  }
  nBndFaces[ib] = ne;

  // Top Edge Faces
  ib = bc2bcList[bcNum[params->create_bcTop]];
  ne = nBndFaces[ib];
  for (int ix=0; ix<nx; ix++) {
    bndFaces[ib][ne][1] = (nx+1)*ny + ix;
    bndFaces[ib][ne][0] = (nx+1)*ny + ix+1;
    bndPts[ib][2*ne]   = bndFaces[ib][ne][0];
    bndPts[ib][2*ne+1] = bndFaces[ib][ne][1];
    ne++;
  }
  nBndFaces[ib] = ne;

  // Left Edge Faces
  ib = bc2bcList[bcNum[params->create_bcLeft]];
  ne = nBndFaces[ib];
  for (int iy=0; iy<ny; iy++) {
    bndFaces[ib][ne][1] = iy*(nx+1);
    bndFaces[ib][ne][0] = (iy+1)*(nx+1);
    bndPts[ib][2*ne]   = bndFaces[ib][ne][0];
    bndPts[ib][2*ne+1] = bndFaces[ib][ne][1];
    ne++;
  }
  nBndFaces[ib] = ne;

  // Right Edge Faces
  ib = bc2bcList[bcNum[params->create_bcRight]];
  ne = nBndFaces[ib];
  for (int iy=0; iy<ny; iy++) {
    bndFaces[ib][ne][0] = iy*(nx+1) + nx;
    bndFaces[ib][ne][1] = (iy+1)*(nx+1) + nx;
    bndPts[ib][2*ne]   = bndFaces[ib][ne][0];
    bndPts[ib][2*ne+1] = bndFaces[ib][ne][1];
    ne++;
  }
  nBndFaces[ib] = ne;

  // Remove duplicates in bndPts
  for (int i=0; i<nBounds; i++) {
    std::sort(bndPts[i], bndPts[i]+bndPts.dim1);
    int* it = std::unique(bndPts[i], bndPts[i]+bndPts.dim1);
    nBndPts[i] = std::distance(bndPts[i],it);
  }
}

void geo::processPeriodicBoundaries(void)
{
  uint nPeriodic, bi, bj, ic;
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
      if (i==j || bndEdges[i]==-1 || bndEdges[j]==-1) continue;
      if (checkPeriodicFaces(e2v[bndEdges[i]],e2v[bndEdges[j]])) {

        /* --- Match found - now take care of transfer from boundary -> internal --- */

        if (i>j) FatalError("How did this happen?!");

        bi = bndEdges[i];
        bj = bndEdges[j];

        // Transfer combined edge from boundary to internal list
        intEdges.push_back(bi);

        // Flag global edge IDs as internal edges
        isBnd[bi] = false;
        isBnd[bj] = false;

        // Fix e2c - add right cell to combined edge, make left cell = -1 in 'deleted' edge
        e2c[bi][1] = e2c[bj][0];
        e2c[bj][0] = -1;

        // Fix c2e - replace 'deleted' edge from right cell with combined edge
        ic = e2c[bi][1];
        int fID = findFirst(c2e[ic],(int)bj,c2ne[ic]);
        c2e[e2c[bi][1]][fID] = bi;

        // Fix c2b - set element-local face to be internal face
        c2b[e2c[bi][1]][fID] = false;

        // Flag edges as gone in boundary edges list
        bndEdges[i] = -1;
        bndEdges[j] = -1;
      }
    }
  }

  // Remove no-longer-existing periodic boundary edges and update nBndEdges
  bndEdges.erase(std::remove(bndEdges.begin(), bndEdges.end(), -1), bndEdges.end());
  nBndEdges = bndEdges.size();
  nFaces = intEdges.size();
}

bool geo::checkPeriodicFaces(int* edge1, int* edge2)
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
