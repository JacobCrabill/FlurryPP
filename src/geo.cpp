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
#include <map>
#include <set>
#include <sstream>
#include <unordered_set>

#include "../include/face.hpp"
#include "../include/intFace.hpp"
#include "../include/boundFace.hpp"
#include "../include/mpiFace.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#include "metis.h"
#endif

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

#ifndef _NO_MPI
  partitionMesh();
#endif

  processConnectivity();

  processPeriodicBoundaries();
}


void geo::processConnectivity()
{
  if (params->rank==0) cout << "Geo: Processing element connectivity" << endl;

  /* --- Setup Edges --- */
  matrix<int> e2v1;
  vector<int> edge(2);

  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2ne[e]; ie++) {  // NOTE: nv may be > ne
      int iep1 = (ie+1)%c2ne[e];
      if (c2v(e,ie) == c2v(e,iep1)) {
        // Collapsed edge - ignore
        continue;
      }
      else if (c2v[e][ie] < c2v(e,iep1)) {
        edge[0] = c2v(e,ie);
        edge[1] = c2v(e,iep1);
      }
      else {
        edge[0] = c2v(e,iep1);
        edge[1] = c2v(e,ie);
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

  /* Flag for whether global face ID corresponds to interior or boundary face
     (note that, at this stage, MPI faces will be considered boundary faces) */
  isBnd.assign(nEdges,0);

  nFaces = 0;
  nBndFaces = 0;
  nMpiFaces = 0;

  for (uint i=0; i<iE.size(); i++) {
    if (iE[i]!=-1) {
      auto ie = findEq(iE,iE[i]);
      if (ie.size()>2) {
        stringstream ss; ss << i;
        string errMsg = "More than 2 cells for edge " + ss.str();
        FatalError(errMsg.c_str());
      }
      else if (ie.size()==2) {
        // Internal Edge which has not yet been added
        intEdges.push_back(iE[i]);
        nFaces++;
      }
      else if (ie.size()==1) {
        // Boundary or MPI Edge
        bndEdges.push_back(iE[i]);
        isBnd[iE[i]] = true;
        nBndFaces++;
      }

      // Mark edges as completed
      vecAssign(iE,ie,-1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bndFaces.resize(nBounds);
  bcType.assign(nBndFaces,-1);
  for (int i=0; i<nBndFaces; i++) {
    int iv1 = e2v(bndEdges[i],0);
    int iv2 = e2v(bndEdges[i],1);
    for (int bnd=0; bnd<nBounds; bnd++) {
      if (findFirst(bndPts[bnd],iv1,bndPts.dim1)!=-1 && findFirst(bndPts[bnd],iv2,bndPts.dim1)!=-1) {
        // The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bndFaces[bnd].insertRow(e2v[bndEdges[i]],INSERT_AT_END,e2v.dim1);
        break;
      }
    }
  }

  /* --- Setup MPI Processor Boundary Faces --- */
#ifndef _NO_MPI
  if (params->nproc > 1) {

    if (params->rank == 0) cout << "Geo: Matching MPI faces" << endl;

    // 1) Get a list of all the MPI faces on the processor
    // These will be all unassigned boundary faces (bcType == -1) - copy over to mpiEdges
    for (int i=0; i<nBndFaces; i++) {
      if (bcType[i] < 0) {
        mpiEdges.push_back(bndEdges[i]);
        bndEdges[i] = -1;
      }
    }
    nMpiFaces = mpiEdges.size();

    // Clean up the bcType and bndEdges arrays now that it's safe to do so [remove mpiFaces from them]
    bndEdges.erase(std::remove(bndEdges.begin(), bndEdges.end(), -1), bndEdges.end());
    bcType.erase(std::remove(bcType.begin(), bcType.end(), -1), bcType.end());
    nBndFaces = bndEdges.size();

    // For future compatibility with 3D mixed meshes: allow faces with different #'s nodes
    // mpi_fptr is like csr matrix ptr (or like eptr from METIS, but for faces instead of eles)
    matrix<int> mpiFaceNodes;
    vector<int> mpiFptr(nMpiFaces+1);
    for (int i=0; i<nMpiFaces; i++) {
      mpiFaceNodes.insertRow(e2v[mpiEdges[i]],INSERT_AT_END,2);
      mpiFptr[i+1] = mpiFptr[i]+2;
    }
    int nMpiNodes = mpiFptr[nMpiFaces];

    // Convert local node ID's to global
    std::transform(mpiFaceNodes.getData(),mpiFaceNodes.getData()+mpiFaceNodes.getSize(),mpiFaceNodes.getData(), [=](int iv){return iv2ivg[iv];} );

    // Get the number of mpiFaces on each processor (for later communication)
    vector<int> nMpiFaces_proc(params->nproc);
    MPI_Allgather(&nMpiFaces,1,MPI_INT,nMpiFaces_proc.data(),1,MPI_INT,MPI_COMM_WORLD);

    // Get the total number of face nodes to be sent on each processor (for later communication)
    vector<int> nMpiNodes_proc(params->nproc);
    MPI_Allgather(&nMpiNodes,1,MPI_INT,nMpiNodes_proc.data(),1,MPI_INT,MPI_COMM_WORLD);

    // 2 for 2D, 4 for 3D; recall that we're treating all elements as being linear, as
    // the extra nodes for quadratic edges or faces are unimportant for determining connectivity
    int maxNodesPerFace = (nDims==2) ? 2 : 4;
    int maxNMpiNodes = getMax(nMpiNodes_proc);
    int maxNMpiFaces = getMax(nMpiFaces_proc);
    matrix<int> mpiFaceNodes_proc(params->nproc,maxNMpiFaces*maxNodesPerFace);
    matrix<int> mpiFptr_proc(params->nproc,maxNMpiFaces+1);
    MPI_Allgather(mpiFaceNodes.getData(),mpiFaceNodes.getSize(),MPI_INT,mpiFaceNodes_proc.getData(),maxNMpiFaces*maxNodesPerFace,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(mpiFptr.data(),mpiFptr.size(),MPI_INT,mpiFptr_proc.getData(),maxNMpiFaces+1,MPI_INT,MPI_COMM_WORLD);

    // Now that we have each processor's boundary nodes, start matching faces
    // Again, note that this is written for to be entirely general instead of 2D-specific
    // Find out what processor each face is adjacent to
    procR.resize(nMpiFaces);
    locF_R.resize(nMpiFaces);
    for (auto &P:procR) P = -1;

    vector<int> tmpFace(maxNodesPerFace);
    vector<int> myFace(maxNodesPerFace);
    for (int p=0; p<params->nproc; p++) {
      if (p == params->rank) continue;

      // Check all of the processor's faces to see if any match our faces
      for (int i=0; i<nMpiFaces_proc[p]; i++) {
        tmpFace.resize(maxNodesPerFace);
        int k = 0;
        for (int j=mpiFptr_proc(p,i); j<mpiFptr_proc(p,i+1); j++) {
          tmpFace[k] = mpiFaceNodes_proc(p,j);
          k++;
        }
        tmpFace.resize(k);

        // See if this face matches any on this processor
        for (int F=0; F<nMpiFaces; F++) {
          if (procR[F] != -1) continue; // Face already matched

          for (int j=0; j<2; j++) {  // crap.  I'll need to set up an 'e2nv' to use here later.
            myFace[j] = mpiFaceNodes(F,j);
          }
          if (compareFaces(myFace,tmpFace)) {
            procR[F] = p;
            locF_R[F] = i;
            break;
          }
        }
      }
    }

    for (auto &P:procR)
      if (P==-1) FatalError("MPI face left unmatched!");

    if (params->rank == 0) cout << "Geo: All MPI faces matched!  nMpiFaces = " << nMpiFaces << endl;
  } // nproc > 1
#endif

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
      if (c2v(ic,j) == c2v(ic,jp1)) {
        // Collapsed edge; ignore
        c2e(ic,j) = -1;
        c2b(ic,j) = 0;
        continue;
      }

      if (c2v(ic,j) < c2v(ic,jp1)) {
        edge[0] = c2v(ic,j);
        edge[1] = c2v(ic,jp1);
      } else {
        edge[0] = c2v(ic,jp1);
        edge[1] = c2v(ic,j);
      }

      auto ie1 = findEq(e2v.getCol(0),edge[0]);
      auto col2 = (e2v.getRows(ie1)).getCol(1);
      int ie2 = findFirst(col2,edge[1]);
      int ie0 = ie1[ie2];

      // Find ID of face within type-specific array
      if (isBnd[ie0]) {
        c2e(ic,j) = ie0;
        c2b(ic,j) = 1;
      }else{
        c2e(ic,j) = ie0;
        c2b(ic,j) = 0;
      }

      if (e2c(ie0,0) == -1) {
        // No cell yet assigned to edge; put on left
        e2c(ie0,0) = ic;
      }else{
        // Put cell on right
        e2c(ie0,1) = ic;
      }
    }
  }
}

void geo::setupElesFaces(vector<ele> &eles, vector<face*> &faces, vector<mpiFace*> &mpiFacesVec)
{
  if (nEles<=0) FatalError("Cannot setup elements array - nEles = 0");

  eles.resize(nEles);
  faces.resize(nFaces+nBndFaces);
  mpiFacesVec.resize(nMpiFaces);

  if (params->rank==0) cout << "Geo: Setting up elements" << endl;

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
      e.nodeID[iv] = c2v(ic,iv);
      e.nodes[iv] = xv[c2v(ic,iv)];
    }

    // Global face IDs for internal & boundary faces
    e.faceID.resize(c2ne[ic]);
    e.bndFace.resize(c2ne[ic]);
    for (int k=0; k<c2ne[ic]; k++) {
      e.bndFace[k] = c2b(ic,k);
      e.faceID[k] = c2e(ic,k);
    }

    ic++;
  }

  /* --- Setup the faces --- */

  vector<int> tmpEdges;

  if (params->rank==0) cout << "Geo: Setting up internal faces" << endl;

  // Internal Faces
  for (int i=0; i<nFaces; i++) {
    intFace *F = new intFace();
    faces[i] = F;
    // Find global face ID of current interior face
    int ie = intEdges[i];
    ic = e2c(ie,0);
    // Find local face ID of global face within first element [on left]
    tmpEdges.assign(c2e[ic],c2e[ic]+c2ne[ic]);
    int fid1 = findFirst(tmpEdges,ie);
    if (e2c(ie,1) == -1) {
      FatalError("Interior face does not have a right element assigned.");
    }else{
      ic = e2c(ie,1);
      tmpEdges.assign(c2e[ic], c2e[ic]+c2ne[ic]);  // List of cell's faces
      int fid2 = findFirst(tmpEdges,ie);           // Which one is this face
      F->initialize(&eles[e2c(ie,0)],&eles[e2c(ie,1)],fid1,fid2,ie,params);
    }
  }

  if (params->rank==0) cout << "Geo: Setting up boundary faces" << endl;

  // Boundary Faces
  for (int i=0; i<nBndFaces; i++) {
    boundFace *B = new boundFace();
    faces[nFaces+i] = B;
    // Find global face ID of current boundary face
    int ie = bndEdges[i];
    ic = e2c(ie,0);
    // Find local face ID of global face within element
    tmpEdges.assign(c2e[ic],c2e[ic]+c2ne[ic]);
    int fid1 = findFirst(tmpEdges,ie);
    if (e2c(ie,1) != -1) {
      FatalError("Boundary face has a right element assigned.");
    }else{
      B->initialize(&eles[e2c(ie,0)],NULL,fid1,bcType[i],ie,params);
    }
  }

#ifndef _NO_MPI
  // MPI Faces
  if (params->nproc > 1) {

    if (params->rank==0) cout << "Geo: Setting up MPI faces" << endl;

    for (int i=0; i<nMpiFaces; i++) {
      mpiFace *F = new mpiFace();
      mpiFacesVec[i] = F;
      // Find global face ID of current boundary face
      int ie = mpiEdges[i];
      ic = e2c(ie,0);
      // Find local face ID of global face within element
      tmpEdges.assign(c2e[ic],c2e[ic]+c2ne[ic]);
      int fid1 = findFirst(tmpEdges,ie);
      if (e2c(ie,1) != -1) {
        FatalError("MPI face has a right element assigned.");
      }else{
        F->procL = params->rank;
        F->procR = procR[i];
        F->initialize(&eles[e2c(ie,0)],NULL,fid1,locF_R[i],i,params);
      }
    }
  }
#endif
}

void geo::readGmsh(string fileName)
{
  ifstream meshFile;
  string str;

  if (params->rank==0) cout << "Geo: Reading mesh file " << fileName << endl;

  meshFile.open(fileName.c_str());
  if (!meshFile.is_open())
    FatalError("Unable to open mesh file.");

  /* --- Read Boundary Conditions & Fluid Field(s) --- */

  // Move cursor to $PhysicalNames
  while(1) {
    getline(meshFile,str);
    if (str.find("$PhysicalNames")!=string::npos) break;
    if(meshFile.eof()) FatalError("$PhysicalNames tag not found in Gmsh file!");
  }

  // Read number of boundaries and fields defined
  int nBnds;              // Temp. variable for # of Gmsh regions ("PhysicalNames")
  meshFile >> nBnds;
  getline(meshFile,str);  // clear rest of line

  nBounds = 0;
  for (int i=0; i<nBnds; i++) {
    string bcStr;
    stringstream ss;
    int bcdim, bcid;

    getline(meshFile,str);
    ss << str;
    ss >> bcdim >> bcid >> bcStr;

    // Remove quotation marks from around boundary condition
    size_t ind = bcStr.find("\"");
    while (ind!=string::npos) {
      bcStr.erase(ind,1);
      ind = bcStr.find("\"");
    }

    // Convert to lowercase to match Flurry's boundary condition strings
    std::transform(bcStr.begin(), bcStr.end(), bcStr.begin(), ::tolower);

    // First, map mesh boundary to boundary condition in input file
    if (!params->meshBounds.count(bcStr)) {
      string errS = "Unrecognized mesh boundary: \"" + bcStr + "\"\n";
      errS += "Boundary names in input file must match those in mesh file.";
      FatalError(errS.c_str());
    }

    bcStr = params->meshBounds[bcStr];

    // Next, check that the requested boundary condition exists
    if (!bcStr2Num.count(bcStr)) {
      string errS = "Unrecognized boundary condition: \"" + bcStr + "\"";
      FatalError(errS.c_str());
    }

    bcList.push_back(bcStr2Num[bcStr]);

    bcIdMap[bcid] = i; // Map Gmsh bcid to Flurry bound index

    if (bcStr.compare("fluid")==0) {
      nDims = bcdim;
      params->nDims = bcdim;
    }
    else {
      nBounds++;
    }
  }

  if (nDims != 2) {
    FatalError("Only 2D meshes are currently supported - check that your mesh is setup properly.");
  }

  /* --- Read Mesh Vertex Locations --- */

  // Move cursor to $Nodes
  meshFile.clear();
  meshFile.seekg(0, ios::beg);
  while(1) {
    getline(meshFile,str);
    if (str.find("$Nodes")!=string::npos) break;
    if(meshFile.eof()) FatalError("$Nodes tag not found in Gmsh file!");
  }

  uint iv;
  meshFile >> nVerts;
  xv.resize(nVerts);
  getline(meshFile,str); // Clear end of line, just in case

  for (int i=0; i<nVerts; i++) {
    meshFile >> iv >> xv[i].x >> xv[i].y >> xv[i].z;
  }

  /* --- Read Element Connectivity --- */

  // Move cursor to $Elements
  meshFile.clear();
  meshFile.seekg(0, ios::beg);
  while(1) {
    getline(meshFile,str);
    if (str.find("$Elements")!=string::npos) break;
    if(meshFile.eof()) FatalError("$Elements tag not found in Gmsh file!");
  }

  int nElesGmsh;
  vector<int> c2v_tmp(9,0);  // Maximum number of nodes/element possible
  vector<set<int>> boundPoints(nBounds);
  map<int,int> eType2nv;
  eType2nv[3] = 4;
  eType2nv[16] = 4;
  eType2nv[10] = 4;

  // Setup bndPts matrix - Just an initial estimate; will be adjusted on the fly
  //bndPts.setup(nBounds,std::round(nNodes/nBounds));
  nBndPts.resize(nBounds);

  // Read total number of interior + boundary elements
  meshFile >> nElesGmsh;
  getline(meshFile,str);    // Clear end of line, just in case

  // For Gmsh node ordering, see: http://geuz.org/gmsh/doc/texinfo/gmsh.html#Node-ordering
  int ic = 0;
  for (int k=0; k<nElesGmsh; k++) {
    int id, eType, nTags, bcid, tmp;
    meshFile >> id >> eType >> nTags;
    meshFile >> bcid;
    bcid = bcIdMap[bcid];
    for (int tag=0; tag<nTags-1; tag++)
      meshFile >> tmp;

    if (bcList[bcid] == NONE) {
      // NOTE: Currently, only quads are supported
      switch(eType) {
      case 2:
        // linear triangle -> linear quad
        c2nv.push_back(4);
        c2ne.push_back(4);
        ctype.push_back(QUAD);
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2];
        c2v_tmp[3] = c2v_tmp[2];
        break;

      case 9:
        // quadratic triangle -> quadratic quad  [corner nodes, then edge-center nodes]
        c2nv.push_back(8);
        c2ne.push_back(4);
        ctype.push_back(QUAD);
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[7];
        c2v_tmp[3] = c2v_tmp[2];
        c2v_tmp[6] = c2v_tmp[2];
        break;

      case 3:
        // linear quadrangle
        c2nv.push_back(4);
        c2ne.push_back(4);
        ctype.push_back(QUAD);
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3];
        break;

      case 16:
        // quadratic 8-node (serendipity) quadrangle
        c2nv.push_back(8);
        c2ne.push_back(4);
        ctype.push_back(QUAD);
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
        break;

      case 10:
        // quadratic (9-node Lagrange) quadrangle (read as 8-node serendipity)
        c2nv.push_back(8);
        c2ne.push_back(4);
        ctype.push_back(QUAD);
        meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7] >> c2v_tmp[8];
        break;

      default:
        cout << "Gmsh Element Type = " << eType << endl;
        FatalError("element type not recognized");
        break;
      }

      // Increase the size of c2v (max # of vertices per cell) if needed
      if (c2v.getDim1()<(uint)c2nv[ic]) {
        for (int dim=c2v.getDim1(); dim<c2nv[ic]; dim++) {
          c2v.addCol();
        }
      }

      // Number of nodes in c2v_tmp may vary, so use pointer rather than vector
      c2v.insertRow(c2v_tmp.data(),-1,c2nv[ic]);

      // Shift every value of c2v by -1 (Gmsh is 1-indexed; we need 0-indexed)
      for(int k=0; k<c2nv[ic]; k++) {
        if(c2v(ic,k)!=0) {
          c2v(ic,k)--;
        }
      }

      ic++;
      getline(meshFile,str); // skip end of line
    }
    else {
      // Boundary cell; put vertices into bndPts
      int nPtsEdge = 0;
      switch(eType) {
      case 1: // Linear edge
        nPtsEdge = 2;
        break;

      case 8: // Quadratic edge
        nPtsEdge = 3;
        break;

      case 26: // Cubic Edge
        nPtsEdge = 4;
        break;

      case 27: // Quartic Edge
        nPtsEdge = 5;
        break;

      case 28: // Quintic Edge
        nPtsEdge = 6;
        break;

      default:
          FatalError("Boundary Element (Face) Type Not Recognized!");
      }

      for (int i=0; i<nPtsEdge; i++) {
        meshFile >> iv;  iv--;
        boundPoints[bcid].insert(iv);  // bcid != bnd index - FIX ME!!!
      }
      getline(meshFile,str);
    }
  } // End of loop over entities

  int maxNBndPts = 0;
  for (int i=0; i<nBounds; i++) {
    nBndPts[i] = boundPoints[i].size();
    maxNBndPts = max(maxNBndPts,nBndPts[i]);
  }

  // Copy temp boundPoints data into bndPts matrix
  bndPts.setup(nBounds,maxNBndPts);
  for (int i=0; i<nBounds; i++) {
    int j = 0;
    for (auto& it:boundPoints[i]) {
      bndPts(i,j) = it;
      j++;
    }
  }

  nEles = c2v.getDim0();

  meshFile.close();
}

void geo::createMesh()
{
  int nx = params->nx;
  int ny = params->ny;
  nDims = params->nDims; // Since I may implement createMesh for 3D later

  if (params->rank==0)
    cout << "Geo: Creating " << nx << "x" << ny << " cartesian mesh" << endl;

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
  bcList.push_back(bcStr2Num[params->create_bcBottom]);
  bcList.push_back(bcStr2Num[params->create_bcRight]);
  bcList.push_back(bcStr2Num[params->create_bcTop]);
  bcList.push_back(bcStr2Num[params->create_bcLeft]);

  // Sort the list & remove any duplicates
  std::sort(bcList.begin(), bcList.end());
  vector<int>::iterator vIt = std::unique(bcList.begin(), bcList.end());
  nBounds = std::distance(bcList.begin(), vIt);     // will I need both an nBounds (i.e., in mesh) and an nBC's (current nBounds)?
  bcList.resize(nBounds);

  // Setup a map so we know where each BC# is inside of bcList
  map<int,int> bc2bcList;

  // Setup boundary connectivity storage
  nFacesPerBnd.assign(nBounds,0);
  bndPts.setup(nBounds,4*nx+4*ny); //(nx+1)*(ny+1));
  nBndPts.resize(nBounds);
  for (int i=0; i<nBounds; i++) {
    bc2bcList[bcList[i]] = i;
  }

  // Bottom Edge Faces
  int ib = bc2bcList[bcStr2Num[params->create_bcBottom]];
  int ne = nFacesPerBnd[ib];
  for (int ix=0; ix<nx; ix++) {
    bndPts[ib][2*ne]   = ix;
    bndPts[ib][2*ne+1] = ix+1;
    ne++;
  }
  nFacesPerBnd[ib] = ne;

  // Top Edge Faces
  ib = bc2bcList[bcStr2Num[params->create_bcTop]];
  ne = nFacesPerBnd[ib];
  for (int ix=0; ix<nx; ix++) {
    bndPts[ib][2*ne]   = (nx+1)*ny + ix+1;
    bndPts[ib][2*ne+1] = (nx+1)*ny + ix;
    ne++;
  }
  nFacesPerBnd[ib] = ne;

  // Left Edge Faces
  ib = bc2bcList[bcStr2Num[params->create_bcLeft]];
  ne = nFacesPerBnd[ib];
  for (int iy=0; iy<ny; iy++) {
    bndPts[ib][2*ne]   = (iy+1)*(nx+1);
    bndPts[ib][2*ne+1] = iy*(nx+1);
    ne++;
  }
  nFacesPerBnd[ib] = ne;

  // Right Edge Faces
  ib = bc2bcList[bcStr2Num[params->create_bcRight]];
  ne = nFacesPerBnd[ib];
  for (int iy=0; iy<ny; iy++) {
    bndPts[ib][2*ne]   = iy*(nx+1) + nx;
    bndPts[ib][2*ne+1] = (iy+1)*(nx+1) + nx;
    ne++;
  }
  nFacesPerBnd[ib] = ne;

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

  for (int i=0; i<nBndFaces; i++) {
    if (bcType[i] == PERIODIC) {
      iPeriodic.push_back(i);
    }
  }

  nPeriodic = iPeriodic.size();

#ifndef _NO_MPI
  if (nPeriodic > 0)
    FatalError("Periodic boundaries not implemented yet with MPI! Recompile for serial.");
#endif

  if (nPeriodic == 0) return;
  if (nPeriodic%2 != 0) FatalError("Expecting even number of periodic faces; have odd number.");
  if (params->rank==0) cout << "Geo: Processing periodic boundaries" << endl;

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
        e2c(bi,1) = e2c[bj][0];
        e2c(bj,0) = -1;

        // Fix c2e - replace 'deleted' edge from right cell with combined edge
        ic = e2c[bi][1];
        int fID = findFirst(c2e[ic],(int)bj,c2ne[ic]);
        c2e(e2c(bi,1),fID) = bi;

        // Fix c2b - set element-local face to be internal face
        c2b(e2c(bi,1),fID) = false;

        // Flag edges as gone in boundary edges list
        bndEdges[i] = -1;
        bndEdges[j] = -1;
      }
    }
  }

  // Remove no-longer-existing periodic boundary edges and update nBndEdges
  bndEdges.erase(std::remove(bndEdges.begin(), bndEdges.end(), -1), bndEdges.end());
  nBndFaces = bndEdges.size();
  nFaces = intEdges.size();
}

bool geo::compareFaces(vector<int> &face1, vector<int> &face2)
{
  uint nv = face1.size();
  if (face2.size() != nv) return false;

  bool found = false;
  // 2D: Check the two possible permuations
  if (nv == 2) {
    if (face1[0] == face2[0] && face1[1] == face2[1]) found = true;
    if (face1[0] == face2[1] && face1[1] == face2[0]) found = true;
  }
  else {
    FatalError("3D not implemented, and expecting linear faces for MPI face-matching.");
  }

  return found;
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


void geo::partitionMesh(void)
{
#ifndef _NO_MPI
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (nproc <= 1) return;

  if (rank == 0) cout << "Geo: Partitioning mesh across " << nproc << " processes" << endl;
  if (rank == 0) cout << "Geo:   Number of elements globally: " << nEles << endl;

  vector<idx_t> eptr(nEles+1);
  vector<idx_t> eind;

  int nn = 0;
  for (int i=0; i<nEles; i++) {
    eind.push_back(c2v(i,0));
    nn++;
    for (int j=1; j<c2nv[i]; j++) {
      if (c2v(i,j)==c2v(i,j-1)) {
        continue; // To deal with collapsed edges
      }
      eind.push_back(c2v(i,j));
      nn++;
    }
    eptr[i+1] = nn;
  }

  int objval;
  vector<int> epart(nEles);
  vector<int> npart(nVerts);

  // int errVal = METIS PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize,
  // idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval,idx_t *epart, idx_t *npart)
  int ncommon = 2; // 2 for 2D, ~3 for 3D

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE; // needed?
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;

  METIS_PartMeshDual(&nEles,&nVerts,eptr.data(),eind.data(),NULL,NULL,
                     &ncommon,&nproc,NULL,options,&objval,epart.data(),npart.data());

  // Copy data to the global arrays & reset local arrays
  nEles_g   = nEles;
  nVerts_g  = nVerts;
  c2v_g     = c2v;      c2v.setup(0,0);
  xv_g      = xv;       xv.resize(0);
  ctype_g   = ctype;    ctype.resize(0);
  c2nv_g    = c2nv;     c2nv.resize(0);
  c2ne_g    = c2ne;     c2ne.resize(0);
  bndPts_g  = bndPts;   bndPts.setup(0,0);
  nBndPts_g = nBndPts;  nBndPts.resize(0);

  // Each processor will now grab its own data according to its rank (proc ID)
  for (int i=0; i<nEles; i++) {
    if (epart[i] == rank) {
      c2v.insertRow(c2v_g[i],-1,c2nv_g[i]);
      ic2icg.push_back(i);
      ctype.push_back(ctype_g[i]);
      c2nv.push_back(c2nv_g[i]);
      c2ne.push_back(c2ne_g[i]);
    }
  }

  nEles = c2v.getDim0();

  // Get list of all vertices (their global IDs) used in new partition
  set<int> myNodes;
  for (int i=0; i<nEles; i++) {
    for (int j=0; j<c2nv[i]; j++) {
      myNodes.insert(c2v(i,j));
    }
  }

  nVerts = myNodes.size();

  // Map from global to local to reset c2v array using local data
  vector<int> ivg2iv(nVerts_g);
  for (auto &iv:ivg2iv) iv = -1;

  // Transfer over all needed vertices to local array
  int nv = 0;
  for (auto &iv: myNodes) {
    xv.push_back(xv_g[iv]);
    iv2ivg.push_back(iv);
    ivg2iv[iv] = nv;
    nv++;
  }

  // bndPts array was already setup globally, so remake keeping only local nodes
  vector<set<int>> boundPoints(nBounds);

  for (int i=0; i<nBounds; i++) {
    for (int j=0; j<nBndPts_g[i]; j++) {
      if (findFirst(iv2ivg,bndPts_g(i,j)) != -1) {
        boundPoints[i].insert(bndPts_g(i,j));
      }
    }
  }

  int maxNBndPts = 0;
  for (int i=0; i<nBounds; i++) {
    nBndPts[i] = boundPoints[i].size();
    maxNBndPts = max(maxNBndPts,nBndPts[i]);
  }

  // Copy temp boundPoints data into bndPts matrix
  // [Transform global node IDs --> local]
  bndPts.setup(nBounds,maxNBndPts);
  for (int i=0; i<nBounds; i++) {
    int j = 0;
    for (auto& it:boundPoints[i]) {
      bndPts(i,j) = ivg2iv[it];
      j++;
    }
  }

  // Lastly, update c2v and bndPts from global --> local node IDs
  std::transform(c2v.getData(),c2v.getData()+c2v.getSize(),c2v.getData(), [=](int ivg){return ivg2iv[ivg];});

  cout << "Geo:   On rank " << rank << ": nEles = " << nEles << endl;

  if (rank == 0) cout << "Geo: Done partitioning mesh" << endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

#include "../include/geo.inl"
