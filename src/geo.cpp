/*!
 * \file geo.cpp
 * \brief Class for handling geometry setup & modification
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

#include "geo.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_set>

#include "face.hpp"
#include "intFace.hpp"
#include "boundFace.hpp"
#include "mpiFace.hpp"
#include "overFace.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#include "metis.h"
#endif

geo::geo()
{
  nodesPerCell = NULL;
}

geo::~geo()
{
  if (nodesPerCell != NULL) {
    delete[] nodesPerCell;
  }
}

void geo::setup(input* params, bool HMG)
{
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;
  meshType = params->meshType;
  gridID = 0;
  gridRank = params->rank;
  nProcGrid = params->nproc;
  rank = params->rank;
  nproc = params->nproc;

  switch(meshType) {
    case READ_MESH:
      readGmsh(params->meshFileName);
      break;

    case CREATE_MESH:
      createMesh();
      break;

#ifndef _NO_MPI
    case OVERSET_MESH:
      // Find out which grid this process will be handling
      nGrids = params->nGrids;
      splitGridProcs();
      readGmsh(params->oversetGrids[gridID]);
      break;
#endif

    default:
      FatalError("Mesh type not recognized.");
  }

  if (HMG)
  {
#ifndef _NO_MPI
    getMpiPartitions();
#endif
    if (nDims == 2)
      processConn2D();
    else
      processConn3D();
  }
  else
  {
#ifndef _NO_MPI
    partitionMesh();
#endif
    processConnectivity();
  }
}

void geo::setup_hmg(input *params, int _gridID, int _gridRank, int _nProcGrid, const vector<int> &_gridIdList, const vector<int> &_epart)
{
  this->params = params;

  nDims = params->nDims;
  nFields = params->nFields;
  meshType = params->meshType;
  rank = params->rank;
  nproc = params->nproc;
  nGrids = params->nGrids;

  gridID = _gridID;
  gridRank = _gridRank;
  nProcGrid = _nProcGrid;
  gridIdList = _gridIdList;

#ifndef _NO_MPI
  if (params->meshType == OVERSET_MESH) {
    // Apperently MPI_Comm objects don't like being copied?  Whatever, just re-create.
    MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);
  } else {
    gridComm = MPI_COMM_WORLD;
  }

  if (nProcGrid > 1 && _epart[0] != -1)
    partitionFromEpart(_epart);
#endif

  processConnectivity();
}

void geo::processConnectivity()
{
  if (params->rank==0) cout << "Geo: Processing element connectivity" << endl;

  if (nDims == 2)
    processConn2D();
  else if (nDims == 3)
    processConn3D();

  processConnExtra();

  processPeriodicBoundaries();

#ifndef _NO_MPI
  /* --- Use TIOGA to find all hole nodes, then setup overset-face connectivity --- */
  if (meshType == OVERSET_MESH) {
    if (nDims == 3)
      setupOverset3D();
    else
      setupOverset2D();

    setFaceIblanks();

    // Since this is initial pre-processing, clear blank/unblanks
    blankCells.clear();
    unblankCells.clear();
  }

  /* --- Setup MPI Processor Boundary Faces --- */
  matchMPIFaces();
#endif

  /* --- Additional setup for moving grids --- */
  if (params->motion) {
    xv0.resize(nVerts);
    for (int i=0; i<nVerts; i++) xv0[i] = point(xv[i],nDims);
    gridVel.setup(nVerts,nDims);
  }
}

void geo::processConn2D(void)
{
  /* --- Setup Edges --- */

  matrix<int> e2v1;
  vector<int> edge(2);

  for (int e=0; e<nEles; e++) {
    for (int ie=0; ie<c2nf[e]; ie++) {  // NOTE: nv may be != ne for 3D
      int iep1 = (ie+1)%c2nf[e];
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
  e2v1.unique(f2v,iE);
  nFaces = f2v.getDim0();

  f2nv.assign(nFaces,2);

  /* --- Generaate Internal and Boundary Face Lists --- */

  /* Flag for whether global face ID corresponds to interior or boundary face
     (note that, at this stage, MPI faces will be considered boundary faces) */
  faceType.assign(nFaces,INTERNAL);

  nIntFaces = 0;
  nBndFaces = 0;
  nMpiFaces = 0;

  intFaces.resize(0);
  bndFaces.resize(0);

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
        intFaces.push_back(iE[i]);
        nIntFaces++;
      }
      else if (ie.size()==1) {
        // Boundary or MPI Edge
        bndFaces.push_back(iE[i]);
        faceType[iE[i]] = BOUNDARY;
        nBndFaces++;
      }

      // Mark edges as completed
      vecAssign(iE,ie,-1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bcFaces.resize(nBounds);
  bcType.assign(nBndFaces,-1);
  for (int i=0; i<nBndFaces; i++) {
    int iv1 = f2v(bndFaces[i],0);
    int iv2 = f2v(bndFaces[i],1);
    for (int bnd=0; bnd<nBounds; bnd++) {
      if (findFirst(bndPts[bnd],iv1,bndPts.dims[1])!=-1 && findFirst(bndPts[bnd],iv2,bndPts.dims[1])!=-1) {
        // The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bcFaces[bnd].insertRow(f2v[bndFaces[i]],INSERT_AT_END,f2v.dims[1]);
        break;
      }
    }
  }

  /* --- Setup Cell-To-Edge, Edge-To-Cell --- */

  c2f.setup(nEles,getMax(c2nf));
  c2b.setup(nEles,getMax(c2nf));
  c2c.setup(nEles,getMax(c2nf));
  c2b.initializeToZero();
  c2c.initializeToValue(-1);
  f2c.setup(nFaces,2);
  f2c.initializeToValue(-1);

  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nf[ic]; j++) {
      int jp1 = (j+1)%(c2nf[ic]);

      // Store edges consistently to allow matching of duplicates
      if (c2v(ic,j) == c2v(ic,jp1)) {
        // Collapsed edge; ignore
        c2f(ic,j) = -1;
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

      auto ie1 = findEq(f2v.getCol(0),edge[0]);
      auto col2 = (f2v.getRows(ie1)).getCol(1);
      int ie2 = findFirst(col2,edge[1]);
      int ie0 = ie1[ie2];

      // Find ID of face within type-specific array
      if (faceType[ie0]>0) {
        c2f(ic,j) = ie0;
        c2b(ic,j) = 1;
      }else{
        c2f(ic,j) = ie0;
        c2b(ic,j) = 0;
      }

      if (f2c(ie0,0) == -1) {
        // No cell yet assigned to edge; put on left
        f2c(ie0,0) = ic;
      }else{
        // Put cell on right
        f2c(ie0,1) = ic;
        // Update c2c for both cells
        int ic2 = f2c(ie0,0);
        vector<int> cellFaces(c2f[ic2],c2f[ic2]+c2nf[ic2]);
        int fid2 = findFirst(cellFaces,ie0);
        c2c(ic,j)     = ic2;
        c2c(ic2,fid2) = ic;
      }
    }
  }
}

void geo::processConn3D(void)
{
  /* --- Setup Single List of All Faces (sorted vertex lists) --- */

  matrix<int> f2v1, e2v1;
  vector<int> f2nv1;

  // Handy map to store local face-vertex lists for each ele type
  map<int,matrix<int>> ct2fv;
  map<int,vector<int>> ct2fnv;
  // --- FIX ORDERING FOR FUTURE USE ---
  ct2fv[HEX].insertRow(vector<int>{0,1,2,3});  // Bottom
  ct2fv[HEX].insertRow(vector<int>{4,5,6,7});  // Top
  ct2fv[HEX].insertRow(vector<int>{3,0,4,7});  // Left
  ct2fv[HEX].insertRow(vector<int>{2,1,5,6});  // Right
  ct2fv[HEX].insertRow(vector<int>{1,0,4,5});  // Front
  ct2fv[HEX].insertRow(vector<int>{3,2,6,7});  // Back
  ct2fnv[HEX] = {4,4,4,4,4,4};
  //ct2fnv[PRISM] = {3,3,4,4,4};
  //ct2fnv[TET] = {3,3,3,3};

  for (int e=0; e<nEles; e++) {
    for (int f=0; f<c2nf[e]; f++) {
      // Get local vertex list for face
      auto iface = ct2fv[ctype[e]].getRow(f);

      // Get global vertex list for face
      vector<int> facev(ct2fnv[ctype[e]][f]);
      for (int i=0; i<ct2fnv[ctype[e]][f]; i++) {
        facev[i] = c2v(e,iface[i]);
        if (i>0 && facev[i] == facev[i-1]) facev[i] = -1;
      }

      // Sort the vertices for easier comparison later
      std::sort(facev.begin(),facev.end());
      f2v1.insertRowUnsized(facev);
      f2nv1.push_back(ct2fnv[ctype[e]][f]);

      e2v1.insertRow({facev[0],facev[1]});
      e2v1.insertRow({facev[1],facev[2]});
      e2v1.insertRow({facev[2],facev[3]});
      e2v1.insertRow({facev[0],facev[3]});
    }
  }

  /* --- Get a unique list of faces --- */

  // NOTE: Could setup f2c here, but I already have another algorithm implemented later

  // iE is of length [original f2v1] with range [final f2v]
  // The number of times a face appears in iF is equal to
  // the number of cells that face touches
  vector<int> iF, iE;
  f2v1.unique(f2v,iF);
  e2v1.unique(e2v,iE);
  nFaces = f2v.getDim0();
  nEdges = e2v.getDim0();

  f2nv.resize(nFaces);
  for (uint i=0; i<f2nv1.size(); i++)
    f2nv[iF[i]] = f2nv1[i];


  /* --- Generate Internal and Boundary Face Lists --- */

  // Flag for whether global face ID corresponds to interior or boundary face
  // (note that, at this stage, MPI faces will be considered boundary faces)
  faceType.assign(nFaces,INTERNAL);

  nIntFaces = 0;
  nBndFaces = 0;
  nMpiFaces = 0;

  for (uint i=0; i<iF.size(); i++) {
    if (iF[i]!=-1) {
      auto ff = findEq(iF,iF[i]);
      if (ff.size()>2) {
        stringstream ss; ss << i;
        string errMsg = "More than 2 cells for face " + ss.str();
        FatalError(errMsg.c_str());
      }
      else if (ff.size()==2) {
        // Internal Edge which has not yet been added
        intFaces.push_back(iF[i]);
        nIntFaces++;
      }
      else if (ff.size()==1) {
        // Boundary or MPI Edge
        bndFaces.push_back(iF[i]);
        faceType[iF[i]] = BOUNDARY;
        nBndFaces++;
      }

      // Mark edges as completed
      vecAssign(iF,ff,-1);
    }
  }

  /* --- Match Boundary Faces to Boundary Conditions --- */

  bcFaces.resize(nBounds);
  bcType.assign(nBndFaces,NONE);
  for (int i=0; i<nBndFaces; i++) {
    for (int bnd=0; bnd<nBounds; bnd++) {
      bool isOnBound = true;
      for (int j=0; j<f2nv[bndFaces[i]]; j++) {
        if (findFirst(bndPts[bnd],f2v(bndFaces[i],j),bndPts.dims[1]) == -1) {
          isOnBound = false;
          break;
        }
      }

      if (isOnBound) {
        //cout << "bndFace matched to bc " << bcList[bnd] << endl;
        // The edge lies on this boundary
        bcType[i] = bcList[bnd];
        bcFaces[bnd].insertRow(f2v[bndFaces[i]],INSERT_AT_END,f2v.dims[1]);
        break;
      }
    }
  }

  /* --- Setup Cell-To-Face, Face-To-Cell --- */

  c2f.setup(nEles,getMax(c2nf));
  c2b.setup(nEles,getMax(c2nf));
  c2c.setup(nEles,getMax(c2nf));
  c2f.initializeToValue(-1);
  c2b.initializeToZero();
  c2c.initializeToValue(-1);
  f2c.setup(nFaces,2);
  f2c.initializeToValue(-1);

  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nf[ic]; j++) {
      // Get local vertex list for face
      auto iface = ct2fv[ctype[ic]].getRow(j);

      // Get global vertex list for face
      int fnv = ct2fnv[ctype[ic]][j];
      vector<int> facev(fnv);
      for (int i=0; i<fnv; i++)
        facev[i] = c2v(ic,iface[i]);

      // Sort the vertices for easier comparison
      std::sort(facev.begin(),facev.end());

      bool found = false;

      // Check if face is actually collapsed (nonexistant) (all nodes identical)
      bool collapsed = true;
      for (int i=1; i<fnv; i++)
        collapsed = ( collapsed && (facev[0]==facev[i]) );

      if (collapsed)
        continue;

      for (int f=0; f<nFaces; f++) {
        if (std::equal(f2v[f],f2v[f]+fnv,facev.begin())) {
          found = true;
          c2f(ic,j) = f;
          break;
        }
      }

      if (!found) FatalError("Unable to match cell face to global face list!");

      int ff = c2f(ic,j);

      // Find ID of face within type-specific array
      if (faceType[ff]>0)
        c2b(ic,j) = 1;
      else
        c2b(ic,j) = 0;

      if (f2c(ff,0) == -1) {
        // No cell yet assigned to edge; put on left
        f2c(ff,0) = ic;
      }else{
        // Put cell on right
        f2c(ff,1) = ic;
        // Update c2c for both cells
        int ic2 = f2c(ff,0);
        vector<int> cellFaces(c2f[ic2],c2f[ic2]+c2nf[ic2]);
        int fid2 = findFirst(cellFaces,ff);
        c2c(ic,j)     = ic2;
        c2c(ic2,fid2) = ic;
      }
    }
  }
}

void geo::processConnExtra(void)
{
  int maxNC;
  if (nDims == 2)
    maxNC = 8;
  else
    maxNC = 26;

  getBoundingBox(xv,minPt,maxPt);

  // Get vertex to vertex/edge connectivity
  vector<set<int>> v2v_tmp(nVerts);
  vector<set<int>> v2e_tmp(nVerts);
  if (nDims == 2) {
    // Faces are the edges
    for (int ie=0; ie<nFaces; ie++) {
      int iv1 = f2v(ie,0);
      int iv2 = f2v(ie,1);
      v2v_tmp[iv1].insert(iv2);
      v2v_tmp[iv2].insert(iv1);
      v2e_tmp[iv1].insert(ie);
      v2e_tmp[iv2].insert(ie);
    }
  } else {
    // Edges were processed separately for 3D
    for (int ie=0; ie<nEdges; ie++) {
      int iv1 = e2v(ie,0);
      int iv2 = e2v(ie,1);
      v2v_tmp[iv1].insert(iv2);
      v2v_tmp[iv2].insert(iv1);
      v2e_tmp[iv1].insert(ie);
      v2e_tmp[iv2].insert(ie);
    }
  }

  v2nv.resize(nVerts);
  for (int iv=0; iv<nVerts; iv++) {
    v2nv[iv] = v2v_tmp[iv].size();
  }

  v2v.setup(nVerts,getMax(v2nv));
  v2e.setup(nVerts,getMax(v2nv));
  for (int iv=0; iv<nVerts; iv++) {
    int j = 0;
    for (auto &iv2:v2v_tmp[iv]) {
      v2v(iv,j) = iv2;
      j++;
    }

    j = 0;
    for (auto &ie:v2e_tmp[iv]) {
      v2e(iv,j) = ie;
      j++;
    }
  }
}

void geo::matchMPIFaces(void)
{
#ifndef _NO_MPI
  if (nProcGrid <= 1) return;

  if (gridRank == 0) {
    if (meshType == OVERSET_MESH)
      cout << "Geo: Grid " << gridID << ": Matching MPI faces" << endl;
    else
      cout << "Geo: Matching MPI faces" << endl;
  }

  // 1) Get a list of all the MPI faces on the processor
  // These will be all unassigned boundary faces (bcType == NONE, or the
  // remaining unmatched bcType == PERIODIC faces)
  // - Copy over to mpiFaces
  mpiPeriodic.resize(0);
  mpiFaces.resize(0);
  for (int i=0; i<nBndFaces; i++) {
    if (bcType[i] <= 0) {
      int ff = bndFaces[i];
      mpiFaces.push_back(ff);
      faceType[ff] = MPI_FACE;
      int periodic = (bcType[i]==PERIODIC) ? 1 : 0;
      mpiPeriodic.push_back(periodic);
      if (nDims == 3) {
        // Get cell ID & cell-local face ID for face-rotation mapping
        mpiCells.push_back(ic2icg[f2c(ff,0)]);
        auto cellFaces = c2f.getRow(f2c(ff,0));
        int fid = findFirst(cellFaces,ff);
        mpiLocF.push_back(fid);
      }
      bndFaces[i] = -1;
    }
  }
  nMpiFaces = mpiFaces.size();

  // Clean up the bcType and bndEdges arrays now that it's safe to do so [remove mpiFaces from them]
  bndFaces.erase(std::remove(bndFaces.begin(), bndFaces.end(), -1), bndFaces.end());
  bcType.erase(std::remove(bcType.begin(), bcType.end(), -1), bcType.end());
  nBndFaces = bndFaces.size();

  // For future compatibility with 3D mixed meshes: allow faces with different #'s nodes
  // mpi_fptr is like csr matrix ptr (or like eptr from METIS, but for faces instead of eles)
  vector<int> mpiFaceNodes;
  vector<int> mpiFptr(nMpiFaces+1);
  for (int i=0; i<nMpiFaces; i++) {
    int ff = mpiFaces[i];
    mpiFaceNodes.insert(mpiFaceNodes.end(),f2v[ff],f2v[ff]+f2nv[ff]);
    mpiFptr[i+1] = mpiFptr[i]+f2nv[ff];
  }

  // Convert local node ID's to global
  std::transform(mpiFaceNodes.begin(),mpiFaceNodes.end(),mpiFaceNodes.begin(), [=](int iv){return iv2ivg[iv];} );

  // Get the number of mpiFaces on each processor (for later communication)
  vector<int> nMpiFaces_proc(nProcGrid);
  MPI_Allgather(&nMpiFaces,1,MPI_INT,nMpiFaces_proc.data(),1,MPI_INT,gridComm);

  // 2 for 2D, 4 for 3D; recall that we're treating all elements as being linear, as
  // the extra nodes for quadratic edges or faces are unimportant for determining connectivity
  int maxNodesPerFace = (nDims==2) ? 2 : 4;
  int maxNMpiFaces = getMax(nMpiFaces_proc);
  matrix<int> mpiFaceNodes_proc(nProcGrid,maxNMpiFaces*maxNodesPerFace);
  matrix<int> mpiFptr_proc(nProcGrid,maxNMpiFaces+1);
  matrix<int> mpiFid_proc(nProcGrid,maxNMpiFaces);

  vector<int> recvCnts(nProcGrid);
  vector<int> recvDisp(nProcGrid);
  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nMpiFaces_proc[i]*maxNodesPerFace;
    recvDisp[i] = i*maxNMpiFaces*maxNodesPerFace;
  }
  MPI_Allgatherv(mpiFaceNodes.data(),mpiFaceNodes.size(),MPI_INT,mpiFaceNodes_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nMpiFaces_proc[i]+1;
    recvDisp[i] = i*(maxNMpiFaces+1);
  }
  MPI_Allgatherv(mpiFptr.data(),mpiFptr.size(),MPI_INT,mpiFptr_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nMpiFaces_proc[i];
    recvDisp[i] = i*maxNMpiFaces;
  }
  MPI_Allgatherv(mpiFaces.data(),nMpiFaces,MPI_INT,mpiFid_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

  matrix<int> mpiCells_proc, mpiLocF_proc;
  if (nDims == 3) {
    // Needed for 3D face-matching (to find relRot)
    mpiCells_proc.setup(nProcGrid,maxNMpiFaces);
    mpiLocF_proc.setup(nProcGrid,maxNMpiFaces);
    MPI_Allgatherv(mpiCells.data(),mpiCells.size(),MPI_INT,mpiCells_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);
    MPI_Allgatherv(mpiLocF.data(),mpiLocF.size(),MPI_INT,mpiLocF_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);
  }

  // For overset meshes, can have an overset face *also* be an MPI face known only to one of the processes
  // So, exchange face Iblank info
  vector<int> mpiIblank;
  vector<int> mpiIblankR;
  matrix<int> mpiIblank_proc;
  if (meshType == OVERSET_MESH) {
    mpiIblank.resize(nMpiFaces);
    mpiIblankR.resize(nMpiFaces);
    mpiIblank_proc.setup(nProcGrid,maxNMpiFaces);
    for (int i=0; i<nMpiFaces; i++)
      mpiIblank[i] = iblankCell[f2c(mpiFaces[i],0)];
      //mpiIblank[i] = iblankFace[mpiFaces[i]];
    MPI_Allgatherv(mpiIblank.data(), nMpiFaces, MPI_INT, mpiIblank_proc.getData(), recvCnts.data(), recvDisp.data(), MPI_INT, gridComm);
  }

  // Now that we have each processor's boundary nodes, start matching faces
  // Again, note that this is written for to be entirely general instead of 2D-specific
  // Find out what processor each face is adjacent to
  procR.resize(nMpiFaces);
  faceID_R.resize(nMpiFaces);
  if (nDims == 3) {
    gIC_R.resize(nMpiFaces);
    mpiLocF_R.resize(nMpiFaces);
  }
  for (auto &P:procR) P = -1;

  vector<int> tmpFace(maxNodesPerFace);
  vector<int> myFace(maxNodesPerFace);
  for (int p=0; p<nProcGrid; p++) {
    if (p == gridRank) continue;

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

        for (int j=0; j<f2nv[mpiFaces[F]]; j++) {
          myFace[j] = mpiFaceNodes[mpiFptr[F]+j];
        }

        bool match;
        if (mpiPeriodic[F])
          match = comparePeriodicMPI(myFace,tmpFace);
        else
          match = compareFaces(myFace,tmpFace);

        if (match) {
          procR[F] = p;
          faceID_R[F] = mpiFid_proc(p,i);
          if (nDims == 3) {
            gIC_R[F] = mpiCells_proc(p,i);
            mpiLocF_R[F] = mpiLocF_proc(p,i);
          }
          if (meshType == OVERSET_MESH)
            mpiIblankR[F] = mpiIblank_proc(p,i);
          break;
        }
      }
    }
  }

  for (auto &P:procR)
    if (P==-1) FatalError("MPI face left unmatched!");

  // For overset grids: Now that we have iblank info for both sides, we can remove
  // any faces which should actually be overset faces
  nMpiFaces = mpiFaces.size();
  if (meshType == OVERSET_MESH) {
    for (int F=0; F<nMpiFaces; F++) {
      if (mpiIblank[F] != NORMAL || mpiIblankR[F] != NORMAL) {
        // Not a normal face; figure out if hole or fringe
        if (mpiIblank[F] == HOLE && mpiIblankR[F] == HOLE)
          iblankFace[mpiFaces[F]] = HOLE;
        else
          iblankFace[mpiFaces[F]] = FRINGE;
      }
    }
  }

  if (params->meshType == OVERSET_MESH) {
    for (auto &ff:mpiFaces)
      for (int j=0; j<f2nv[ff]; j++)
        mpiNodes.insert(f2v(ff,j));
  }

  int nFacesTotalGrid;
  MPI_Allreduce(&nMpiFaces,&nFacesTotalGrid,1,MPI_INT,MPI_SUM,gridComm);

  if (gridRank == 0)
    cout << "Geo: Grid " << gridID << ": All MPI faces matched!  nMpiFaces = " << nFacesTotalGrid/2 << endl;
#endif
}

void geo::setupElesFaces(input *params, vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mpiFacesVec, vector<shared_ptr<overFace>> &overFacesVec)
{
  if (nEles<=0) FatalError("Cannot setup elements array - nEles = 0");

  if (gridRank==0) {
    if (meshType == OVERSET_MESH)
      cout << "Geo: Grid " << gridID << ": Setting up elements" << endl;
    else
      cout << "Geo: Setting up elements" << endl;
  }

  eles.resize(0);
  faces.resize(0);
  mpiFacesVec.resize(0);
  overFacesVec.resize(0);

  // Setup the elements

  eleMap.assign(nEles,-1);

  int nc = 0;
  for (int ic=0; ic<nEles; ic++) {
    // Skip any hole cells
    if (meshType == OVERSET_MESH && iblankCell[ic] == HOLE) continue;

    shared_ptr<ele> e = make_shared<ele>();
    e->ID = ic;
    if (nProcGrid>1)
      e->IDg = ic2icg[ic];
    else
      e->IDg = ic;
    e->eType = ctype[ic];
    e->nNodes = c2nv[ic];
    if (nDims == 2)
      e->nMpts = 4;
    else
      e->nMpts = 8;

    // Shape [mesh] nodes
    e->nodeID.resize(c2nv[ic]);
    e->nodes.resize(c2nv[ic]);
    for (int iv=0; iv<c2nv[ic]; iv++) {
      e->nodeID[iv] = c2v(ic,iv);
      e->nodes[iv] = point(xv[c2v(ic,iv)],nDims);
    }

    // Global face IDs for internal & boundary faces
    e->faceID.resize(c2nf[ic]);
    e->bndFace.resize(c2nf[ic]);
    for (int k=0; k<c2nf[ic]; k++) {
      e->bndFace[k] = c2b(ic,k);
      e->faceID[k] = c2f(ic,k);
    }

    eles.push_back(e);
    eleMap[ic] = nc;
    nc++;
  }

  /* --- Setup the faces --- */

  if (meshType == OVERSET_MESH ) {
    /* --- Get a unique, sorted list of all overset faces (either from Gmsh
     * boundary condition, or from TIOGA-based cell blanking) --- */

    for (int ff=0; ff<nFaces; ff++) if (iblankFace[ff] == FRINGE) overFaces.insert(ff);
    for (int ff=0; ff<nBndFaces; ff++) if (bcType[ff] == OVERSET) overFaces.insert(bndFaces[ff]);
  }

  vector<int> cellFaces;

  if (gridRank==0) {
    if (meshType == OVERSET_MESH)
      cout << "Geo: Grid " << gridID << ": Setting up internal faces" << endl;
    else
      cout << "Geo: Setting up internal faces" << endl;
  }

  faceMap.assign(nFaces,-1);
  currFaceType.assign(nFaces,HOLE_FACE);

  // Internal Faces
  for (auto &ff: intFaces) {
    // Skip any hole faces
    if (meshType == OVERSET_MESH && iblankFace[ff] != NORMAL) continue;

    shared_ptr<face> iface = make_shared<intFace>();

    int ic1 = f2c(ff,0);
    // Find local face ID of global face within first element [on left]
    cellFaces.assign(c2f[ic1],c2f[ic1]+c2nf[ic1]);
    int fid1 = findFirst(cellFaces,ff);
    if (f2c(ff,1) == -1) {
      FatalError("Interior face does not have a right element assigned.");
    }
    else {
      int ic2 = f2c(ff,1);
      cellFaces.assign(c2f[ic2], c2f[ic2]+c2nf[ic2]);  // List of cell's faces
      int fid2 = findFirst(cellFaces,ff);           // Which one is this face
      int relRot = compareOrientation(ic1,fid1,ic2,fid2);
      struct faceInfo info;
      info.IDR = fid2;
      info.relRot = relRot;
      ic1 = eleMap[ic1];
      ic2 = eleMap[ic2];
      if (ic2==-1) FatalError("Internal face has right cell blanked.");
      iface->initialize(eles[ic1],eles[ic2],ff,fid1,info,params);
    }

    faces.push_back(iface);
    faceMap[ff] = faces.size()-1;
    currFaceType[ff] = INTERNAL;
  }
  // New # of internal faces (after excluding blanked faces)
  nIntFaces = faces.size();

  if (gridRank==0) {
    if (meshType == OVERSET_MESH)
      cout << "Geo: Grid " << gridID << ": Setting up boundary faces" << endl;
    else
      cout << "Geo: Setting up boundary faces" << endl;
  }

  // Boundary Faces
  for (int i=0; i<nBndFaces; i++) {
    // Find global face ID of current boundary face
    int ff = bndFaces[i];

    if ( (meshType == OVERSET_MESH && iblankFace[ff] != NORMAL) || bcType[i] == OVERSET)
      continue;

    shared_ptr<face> bface = make_shared<boundFace>();

    int ic = f2c(ff,0);
    // Find local face ID of global face within element
    cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
    int fid1 = findFirst(cellFaces,ff);
    if (f2c(ff,1) != -1) {
      FatalError("Boundary face has a right element assigned.");
    }else{
      struct faceInfo info;
      info.bcType = bcType[i];
      info.isBnd = 1;
      ic = eleMap[ic];
      shared_ptr<ele> nullEle;  // Since just giving the funciton 'NULL' isn't possible
      bface->initialize(eles[ic],nullEle,ff,fid1,info,params);
    }

    faces.push_back(bface);
    faceMap[ff] = faces.size()-1;
    currFaceType[ff] = BOUNDARY;
  }
  // New # of boundary faces (after excluding blanked faces)
  nBndFaces = faces.size() - nIntFaces;

#ifndef _NO_MPI
  // MPI Faces
  if (params->nproc > 1) {

    if (gridRank==0) {
      if (meshType == OVERSET_MESH)
        cout << "Geo: Grid " << gridID << ": Setting up MPI faces" << endl;
      else
        cout << "Geo: Setting up MPI faces" << endl;
    }

    for (int i=0; i<nMpiFaces; i++) {
      // Find global face ID of current boundary face
      int ff = mpiFaces[i];

      if (meshType == OVERSET_MESH && iblankFace[ff] != NORMAL) continue;

      shared_ptr<mpiFace> mface = make_shared<mpiFace>();

      int ic = f2c(ff,0);
      // Find local face ID of global face within element
      int fid1;
      if (nDims == 2) {
        cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
        fid1 = findFirst(cellFaces,ff);
      }
      else {
        fid1 = mpiLocF[i];
      }
      if (f2c(ff,1) != -1) {
        FatalError("MPI face has a right element assigned.");
      }else{
        int relRot = 0;
        if (nDims == 3) {
          // Find the relative orientation (rotation) between left & right faces
          relRot = compareOrientationMPI(ic,fid1,gIC_R[i],mpiLocF_R[i],mpiPeriodic[i]);
        }
        struct faceInfo info;
        info.IDR = faceID_R[i];
        info.relRot = relRot;
        info.procL = gridRank;
        info.procR = procR[i];
        info.isMPI = 1;
        info.gridComm = gridComm;  // Note that this is equivalent to MPI_COMM_WORLD if non-overset (ngrids = 1)
        ic = eleMap[ic];
        shared_ptr<ele> nullEle;  // Since just giving the funciton 'NULL' isn't possible
        mface->initialize(eles[ic],nullEle,ff,fid1,info,params);
      }

      mpiFacesVec.push_back(mface);
      faceMap[ff] = mpiFacesVec.size()-1;
      currFaceType[ff] = MPI_FACE;
    }
  }
#endif

  // Overset Faces
  if (meshType == OVERSET_MESH) {
    set<int> rm_ids;
    for (auto &ff: overFaces) {
      iblankFace[ff] = FRINGE; // To ensure consistency

      shared_ptr<overFace> oface = make_shared<overFace>();

      int ic = f2c(ff,0);
      if (ic == -1 || iblankCell[f2c(ff,0)] == HOLE) {
        if (f2c(ff,1) == -1 || iblankCell[f2c(ff,1)] == HOLE) {
          // This happens when a fringe face is ALSO an MPI-boundary face
          // Since the other processor has the non-blanked cell, just ignore the face here
          iblankFace[ff] = HOLE;
          rm_ids.insert(ff);
          continue;
        }
        ic = f2c(ff,1);
      }

      // Find local face ID of global face within first element [on left]
      cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
      int fid = findFirst(cellFaces,ff);

      struct faceInfo info;
      ic = eleMap[ic];
      shared_ptr<ele> nullEle;  // Since just giving the funciton 'NULL' isn't possible
      oface->initialize(eles[ic],nullEle,ff,fid,info,params);

      overFacesVec.push_back(oface);
      faceMap[ff] = overFacesVec.size()-1;
      currFaceType[ff] = OVER_FACE;
    }

    for (auto &ff:rm_ids) overFaces.erase(ff);
  }
}

void geo::readGmsh(string fileName)
{
  ifstream meshFile;
  string str;

  if (meshType == OVERSET_MESH) {
    if (gridRank==0) cout << "Geo: Grid " << gridID << ": Reading mesh file " << fileName << endl;
  }
  else {
    if (gridRank==0) cout << "Geo: Reading mesh file " << fileName << endl;
  }

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

    // Map the Gmsh PhysicalName to the input-file-specified Flurry boundary condition
    bcStr = params->meshBounds[bcStr];

    // Next, check that the requested boundary condition exists
    if (!bcStr2Num.count(bcStr)) {
      string errS = "Unrecognized boundary condition: \"" + bcStr + "\"";
      FatalError(errS.c_str());
    }

    if (bcStr.compare("fluid")==0) {
      nDims = bcdim;
      params->nDims = bcdim;
      bcIdMap[bcid] = -1;
    }
    else {
      bcList.push_back(bcStr2Num[bcStr]);
      bcIdMap[bcid] = nBounds; // Map Gmsh bcid to Flurry bound index
      nBounds++;
    }
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
  xv.setup(nVerts,nDims);
  getline(meshFile,str); // Clear end of line, just in case

  for (int i=0; i<nVerts; i++) {
    meshFile >> iv >> xv(i,0) >> xv(i,1);
    if (nDims == 3) meshFile >> xv(i,2);
    getline(meshFile,str);
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
  vector<int> c2v_tmp(27,0);  // Maximum number of nodes/element possible
  vector<set<int>> boundPoints(nBounds);
  map<int,int> eType2nv;
  eType2nv[3] = 4;  // Linear quad
  eType2nv[16] = 4; // Quadratic serendipity quad
  eType2nv[10] = 4; // Quadratic Lagrange quad
  eType2nv[8] = 8;  // Linear hex

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

    if (bcid == -1) {
      // NOTE: Currently, only quads are supported
      switch(eType) {
        case 2:
          // linear triangle -> linear quad
          c2nv.push_back(4);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2];
          c2v_tmp[3] = c2v_tmp[2];
          break;

        case 9:
          // quadratic triangle -> quadratic quad  [corner nodes, then edge-center nodes]
          c2nv.push_back(8);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[7];
          c2v_tmp[3] = c2v_tmp[2];
          c2v_tmp[6] = c2v_tmp[2];
          break;

        case 3:
          // linear quadrangle
          c2nv.push_back(4);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3];
          break;

        case 16:
          // quadratic 8-node (serendipity) quadrangle
          c2nv.push_back(8);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
          break;

        case 10:
          // quadratic (9-node Lagrange) quadrangle (read as 8-node serendipity)
          c2nv.push_back(8);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
          break;

        case 36:
          // cubic (16-node Lagrange) quadrangle
          c2nv.push_back(16);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 16; i++) meshFile >> c2v_tmp[i];
          break;

        case 37:
          // quartic (25-node Lagrange) quadrangle
          c2nv.push_back(25);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 25; i++) meshFile >> c2v_tmp[i];
          break;

        case 38:
          // quintic (36-node Lagrange) quadrangle
          c2nv.push_back(36);
          c2v_tmp.resize(36);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 36; i++) meshFile >> c2v_tmp[i];
          break;

        case 47:
          // 6th-order 49-node Lagrange quadrangle
          c2nv.push_back(49);
          c2v_tmp.resize(49);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 49; i++) meshFile >> c2v_tmp[i];
          break;

        case 48:
          // 7th-order 64-node Lagrange quadrangle
          c2nv.push_back(64);
          c2v_tmp.resize(64);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 64; i++) meshFile >> c2v_tmp[i];
          break;

        case 49:
          // 8th-order 81-node Lagrange quadrangle
          c2nv.push_back(81);
          c2v_tmp.resize(81);
          c2nf.push_back(4);
          ctype.push_back(QUAD);
          for (int i = 0; i < 81; i++) meshFile >> c2v_tmp[i];
          break;

        case 5:
          // Linear hexahedron
          c2nv.push_back(8);
          c2nf.push_back(6);
          ctype.push_back(HEX);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
          break;

        case 17:
          // Quadratic (20-Node Serendipity) Hexahedron
          c2nv.push_back(20);
          c2nf.push_back(6);
          ctype.push_back(HEX);
          // Corner Nodes
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
          // Edge Nodes
          meshFile >> c2v_tmp[8] >> c2v_tmp[11] >> c2v_tmp[12] >> c2v_tmp[9] >> c2v_tmp[13] >> c2v_tmp[10];
          meshFile >> c2v_tmp[14] >> c2v_tmp[15] >> c2v_tmp[16] >> c2v_tmp[19] >> c2v_tmp[17] >> c2v_tmp[18];
          break;

        case 12:
          // Quadratic (27-Node Lagrange) Hexahedron (read as 20-node serendipity)
          c2nv.push_back(20);
          c2nf.push_back(6);
          ctype.push_back(HEX);
          c2v_tmp.resize(20);
          // Corner Nodes
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[3] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6] >> c2v_tmp[7];
          // Edge Nodes
          meshFile >> c2v_tmp[8] >> c2v_tmp[11] >> c2v_tmp[12] >> c2v_tmp[9] >> c2v_tmp[13] >> c2v_tmp[10];
          meshFile >> c2v_tmp[14] >> c2v_tmp[15] >> c2v_tmp[16] >> c2v_tmp[19] >> c2v_tmp[17] >> c2v_tmp[18];
          break;

        case 4:
          // Linear tetrahedron; read as collapsed-face hex
          c2nv.push_back(4);
          c2nf.push_back(4);
          ctype.push_back(HEX);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[4];
          c2v_tmp[3] = 2;
          c2v_tmp[5] = c2v_tmp[4];
          c2v_tmp[6] = c2v_tmp[4];
          c2v_tmp[6] = c2v_tmp[4];
          break;

        case 6:
          // Linear prism; read as collapsed-face hex
          c2nv.push_back(8);
          c2nf.push_back(6);
          ctype.push_back(HEX);
          meshFile >> c2v_tmp[0] >> c2v_tmp[1] >> c2v_tmp[2] >> c2v_tmp[4] >> c2v_tmp[5] >> c2v_tmp[6];
          c2v_tmp[3] = c2v_tmp[2];
          c2v_tmp[7] = c2v_tmp[6];
          break;

        default:
          cout << "Gmsh element ID " << k << ", Gmsh Element Type = " << eType << endl;
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
      int nPtsFace = 0;
      switch(eType) {
        case 1: // Linear edge
          nPtsFace = 2;
          break;

        case 2: // Linear triangle
          nPtsFace = 3;
          break;

        case 3: // Linear quad
          nPtsFace = 4;
          break;

        case 10: // Quadratic (Lagrange) quad
          nPtsFace = 4;
          break;

        case 16: // Quadratic (Serendipity) quad
          nPtsFace = 4;
          break;

        case 8: // Quadratic edge
          nPtsFace = 3;
          break;

        case 26: // Cubic Edge
          nPtsFace = 4;
          break;

        case 27: // Quartic Edge
          nPtsFace = 5;
          break;

        case 28: // Quintic Edge
          nPtsFace = 6;
          break;

        case 64: // Order 8
          nPtsFace = 9;
          break;

        default:
          cout << "Gmsh element ID " << k << ", Gmsh Element Type = " << eType << endl;
          FatalError("Boundary Element (Face) Type Not Recognized!");
      }

      for (int i=0; i<nPtsFace; i++) {
        meshFile >> iv;  iv--;
        boundPoints[bcid].insert(iv);
      }
      getline(meshFile,str);
    }
  } // End of loop over entities

  nNodesPerCell = getMax(c2nv);

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
  int nz = params->nz;
  nDims = params->nDims;

  if (nDims == 2)
    nz = 1;

  if (params->rank==0)
    cout << "Geo: Creating " << nx << "x" << ny << "x" << nz << " cartesian mesh" << endl;

  double xmin = params->xmin;
  double xmax = params->xmax;
  double ymin = params->ymin;
  double ymax = params->ymax;
  double zmin = params->zmin;
  double zmax = params->zmax;

  double dx = (xmax-xmin)/nx;
  double dy = (ymax-ymin)/ny;
  double dz = (zmax-zmin)/nz;

  params->periodicDX = xmax-xmin;
  params->periodicDY = ymax-ymin;
  params->periodicDZ = zmax-zmin;

  nEles = nx*ny*nz;
  vector<int> c2v_tmp;

  if (nDims == 2) {
    nVerts = (nx+1)*(ny+1);
    xv.setup(nVerts,nDims);

    c2nv.assign(nEles,4);
    c2nf.assign(nEles,4);
    ctype.assign(nEles,QUAD);

    /* --- Setup Vertices --- */

    c2v_tmp.assign(4,0);

    int nv = 0;

    for (int i=0; i<ny+1; i++) {
      for (int j=0; j<nx+1; j++) {
        xv(nv,0) = xmin + j*dx;
        xv(nv,1) = ymin + i*dy;
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
  }
  else if (nDims == 3) {
    nVerts = (nx+1)*(ny+1)*(nz+1);
    xv.setup(nVerts,nDims);

    c2nv.assign(nEles,8);
    c2nf.assign(nEles,6);
    ctype.assign(nEles,HEX);

    /* --- Setup Vertices --- */

    c2v_tmp.assign(8,0);

    int nv = 0;

    for (int k=0; k<nz+1; k++) {
      for (int j=0; j<ny+1; j++) {
        for (int i=0; i<nx+1; i++) {
          xv(nv,0) = xmin + i*dx;
          xv(nv,1) = ymin + j*dy;
          xv(nv,2) = zmin + k*dz;
          nv++;
        }
      }
    }

    /* --- Setup Elements --- */

    for (int k=0; k<nz; k++) {
      for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
          c2v_tmp[0] = i + (nx+1)*(j   + (ny+1)*k);
          c2v_tmp[1] = i + (nx+1)*(j   + (ny+1)*k) + 1;
          c2v_tmp[2] = i + (nx+1)*(j+1 + (ny+1)*k) + 1;
          c2v_tmp[3] = i + (nx+1)*(j+1 + (ny+1)*k);

          c2v_tmp[4] = i + (nx+1)*(j   + (ny+1)*(k+1));
          c2v_tmp[5] = i + (nx+1)*(j   + (ny+1)*(k+1)) + 1;
          c2v_tmp[6] = i + (nx+1)*(j+1 + (ny+1)*(k+1)) + 1;
          c2v_tmp[7] = i + (nx+1)*(j+1 + (ny+1)*(k+1));
          c2v.insertRow(c2v_tmp);
        }
      }
    }
  }

  /* --- Setup Boundaries --- */

  // List of all boundary conditions being used (bcNum maps string->int)
  bcList.push_back(bcStr2Num[params->create_bcBottom]);
  bcList.push_back(bcStr2Num[params->create_bcRight]);
  bcList.push_back(bcStr2Num[params->create_bcTop]);
  bcList.push_back(bcStr2Num[params->create_bcLeft]);

  if (nDims == 3) {
    bcList.push_back(bcStr2Num[params->create_bcFront]);
    bcList.push_back(bcStr2Num[params->create_bcBack]);
  }

  // Sort the list & remove any duplicates
  std::sort(bcList.begin(), bcList.end());
  vector<int>::iterator vIt = std::unique(bcList.begin(), bcList.end());
  nBounds = std::distance(bcList.begin(), vIt);     // will I need both an nBounds (i.e., in mesh) and an nBC's (current nBounds)?
  bcList.resize(nBounds);

  // Setup a map so we know where each BC# is inside of bcList
  map<int,int> bc2bcList;

  // Setup boundary connectivity storage
  nFacesPerBnd.assign(nBounds,0);
  if (nDims == 2) {
    bndPts.setup(nBounds,2*4*(std::max(nx,ny)+1));
  }
  else if (nDims == 3) {
    int maxN_BFace = std::max(nx*ny,nx*nz);
    maxN_BFace = std::max(maxN_BFace,ny*nz);
    bndPts.setup(nBounds,6*4*maxN_BFace);
  }
  nBndPts.resize(nBounds);
  for (int i=0; i<nBounds; i++) {
    bc2bcList[bcList[i]] = i;
  }

  /* --- Setup Boundary Faces --- */

  if (nDims == 2) {
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
  }
  else if (nDims == 3) {
    // Bottom Side Faces
    int ib = bc2bcList[bcStr2Num[params->create_bcBottom]];
    int nf = nFacesPerBnd[ib];
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
        bndPts(ib,4*nf)   = iy*(nx+1) + ix;
        bndPts(ib,4*nf+1) = iy*(nx+1) + ix + 1;
        bndPts(ib,4*nf+2) = (iy+1)*(nx+1) + ix + 1;
        bndPts(ib,4*nf+3) = (iy+1)*(nx+1) + ix;
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;

    // Top Side Faces
    ib = bc2bcList[bcStr2Num[params->create_bcTop]];
    nf = nFacesPerBnd[ib];
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
        bndPts(ib,4*nf)   = (nx+1)*(ny+1)*nz + iy*(nx+1) + ix;
        bndPts(ib,4*nf+1) = (nx+1)*(ny+1)*nz + iy*(nx+1) + ix+1;
        bndPts(ib,4*nf+2) = (nx+1)*(ny+1)*nz + (iy+1)*(nx+1) + ix+1;
        bndPts(ib,4*nf+3) = (nx+1)*(ny+1)*nz + (iy+1)*(nx+1) + ix;
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;

    // Left Side Faces (x = xmin)
    ib = bc2bcList[bcStr2Num[params->create_bcLeft]];
    nf = nFacesPerBnd[ib];
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        bndPts(ib,4*nf)   = iz*(nx+1)*(ny+1) + iy*(nx+1);
        bndPts(ib,4*nf+1) = (iz+1)*(nx+1)*(ny+1) + iy*(nx+1);
        bndPts(ib,4*nf+2) = (iz+1)*(nx+1)*(ny+1) + (iy+1)*(nx+1);
        bndPts(ib,4*nf+3) = iz*(nx+1)*(ny+1) + (iy+1)*(nx+1);
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;

    // Right Side Faces (x = xmax)
    ib = bc2bcList[bcStr2Num[params->create_bcRight]];
    nf = nFacesPerBnd[ib];
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        bndPts(ib,4*nf)   = iz*(nx+1)*(ny+1) + iy*(nx+1) + nx;
        bndPts(ib,4*nf+1) = (iz+1)*(nx+1)*(ny+1) + iy*(nx+1) + nx;
        bndPts(ib,4*nf+2) = (iz+1)*(nx+1)*(ny+1) + (iy+1)*(nx+1) + nx;
        bndPts(ib,4*nf+3) = iz*(nx+1)*(ny+1) + (iy+1)*(nx+1) + nx;
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;


    // Back Side Faces (y = ymin)
    ib = bc2bcList[bcStr2Num[params->create_bcBack]];
    nf = nFacesPerBnd[ib];
    for (int iz=0; iz<nz; iz++) {
      for (int ix=0; ix<nx; ix++) {
        bndPts(ib,4*nf)   = iz*(nx+1)*(ny+1) + ix;
        bndPts(ib,4*nf+1) = (iz+1)*(nx+1)*(ny+1) + ix;
        bndPts(ib,4*nf+2) = (iz+1)*(nx+1)*(ny+1) + ix + 1;
        bndPts(ib,4*nf+3) = iz*(nx+1)*(ny+1) + ix + 1;
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;

    // Front Side Faces (y = ymax)
    ib = bc2bcList[bcStr2Num[params->create_bcFront]];
    nf = nFacesPerBnd[ib];
    for (int iz=0; iz<nz; iz++) {
      for (int ix=0; ix<nx; ix++) {
        bndPts(ib,4*nf)   = iz*(nx+1)*(ny+1) + ny*(nx+1) + ix;
        bndPts(ib,4*nf+1) = (iz+1)*(nx+1)*(ny+1) + ny*(nx+1) + ix;
        bndPts(ib,4*nf+2) = (iz+1)*(nx+1)*(ny+1) + ny*(nx+1) + ix + 1;
        bndPts(ib,4*nf+3) = iz*(nx+1)*(ny+1) + ny*(nx+1) + ix + 1;
        nf++;
      }
    }
    nFacesPerBnd[ib] = nf;
  }

  // Remove duplicates in bndPts
  for (int i=0; i<nBounds; i++) {
    std::sort(bndPts[i], bndPts[i]+bndPts.dims[1]);
    int* it = std::unique(bndPts[i], bndPts[i]+bndPts.dims[1]);
    nBndPts[i] = std::distance(bndPts[i],it);
  }
  int maxNBndPts = getMax(nBndPts);
  bndPts.removeCols(bndPts.dims[1]-maxNBndPts);
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

  if (nPeriodic == 0) return;
  if (nPeriodic%2 != 0 && nProcGrid==1) FatalError("Expecting even number of periodic faces; have odd number.");
  if (params->rank==0) cout << "Geo: Processing periodic boundaries" << endl;

  int nUnmatched = 0;

  for (auto& i:iPeriodic) {
    if (bndFaces[i]==-10) continue;
    bool match = false;
    for (auto& j:iPeriodic) {
      if (i==j || bndFaces[i]==-10 || bndFaces[j]==-10) continue;
      if (nDims == 2) {
        match = checkPeriodicFaces(f2v[bndFaces[i]],f2v[bndFaces[j]]);
      }
      else {
        auto face1 = f2v.getRow(bndFaces[i]);
        auto face2 = f2v.getRow(bndFaces[j]);
        match = checkPeriodicFaces3D(face1, face2);
      }
      if (match) {

        /* --- Match found - now take care of transfer from boundary -> internal --- */

        if (i>j) FatalError("How did this happen?!");

        bi = bndFaces[i];
        bj = bndFaces[j];

        // Transfer combined edge from boundary to internal list
        intFaces.push_back(bi);

        // Flag global edge IDs as internal faces
        faceType[bi] = INTERNAL;
        faceType[bj] = INTERNAL;

        // Fix f2c - add right cell to combined face, make left cell = -1 in 'deleted' face
        f2c(bi,1) = f2c(bj,0);
        f2c(bj,0) = -1;

        // Fix c2f - replace 'deleted' edge from right cell with combined face
        ic = f2c[bi][1];
        int fID = findFirst(c2f[ic],(int)bj,c2nf[ic]);
        c2f(f2c(bi,1),fID) = bi;

        // Fix c2b - set element-local face to be internal face
        c2b(f2c(bi,1),fID) = false;

        // Flag edges as gone in boundary edges list
        bndFaces[i] = -10;
        bndFaces[j] = -10;
        bcType[i] = -10;
        bcType[j] = -10;

        break;
      }
    }

    if (!match)
      nUnmatched++;
  }

  // Remove no-longer-existing periodic boundary faces and update nBndFaces
  bndFaces.erase(std::remove(bndFaces.begin(), bndFaces.end(), -10), bndFaces.end());
  bcType.erase(std::remove(bcType.begin(), bcType.end(), -10), bcType.end());
  nBndFaces = bndFaces.size();
  nIntFaces = intFaces.size();

#ifdef _NO_MPI
  // Can't remember why, but this doesn't apply for MPI cases
  if (nUnmatched>0)
    FatalError("Unmatched periodic faces exist.");
#endif
}

bool geo::compareFaces(vector<int> &face1, vector<int> &face2)
{
  uint nv = face1.size();
  if (face2.size() != nv) return false;

  std::sort(face1.begin(),face1.end());
  std::sort(face2.begin(),face2.end());

  bool found = true;
  for (uint i=0; i<nv; i++) {
    if (face1[i] != face2[i]) found = false;
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

  // None of the above
  return false;
}

bool geo::checkPeriodicFaces3D(vector<int> &face1, vector<int> &face2)
{
  if (face1.size() != face2.size())
    return false;

  double tol = params->periodicTol;
  double dx = params->periodicDX;
  double dy = params->periodicDY;
  double dz = params->periodicDZ;

  /* --- Compare faces using normal vectors: normals should be aligned
   * and offset by norm .dot. {dx,dy,dz} --- */

  Vec3 vec1, vec2;

  // Calculate face normal & centriod for face 1
  Vec3 norm1;
  point c1;
  vec1 = point(xv[face1[1]]) - point(xv[face1[0]]);
  vec2 = point(xv[face1[2]]) - point(xv[face1[0]]);
  for (uint j=0; j<face1.size(); j++)
    c1 += xv[face1[j]];
  c1 /= face1.size();

  norm1[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  norm1[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  norm1[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  // Normalize
  double magNorm1 = sqrt(norm1[0]*norm1[0]+norm1[1]*norm1[1]+norm1[2]*norm1[2]);
  norm1 /= magNorm1;

  // Calculate face normal & centroid for face 2
  Vec3 norm2;
  point c2;
  vec1 = point(xv[face2[1]]) - point(xv[face2[0]]);
  vec2 = point(xv[face2[2]]) - point(xv[face2[0]]);
  for (uint j=0; j<face2.size(); j++)
    c2 += point(xv[face2[j]]);
  c2 /= face2.size();
  norm2[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  norm2[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  norm2[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  // Normalize
  double magNorm2 = sqrt(norm2[0]*norm2[0]+norm2[1]*norm2[1]+norm2[2]*norm2[2]);
  norm2 /= magNorm2;

  // Check for same orientation - norm1 .dot. norm2 should be +/- 1
  double dot = norm1*norm2;
  if (abs(1-abs(dot))>tol) return false;

  // Check offset distance - norm .times. {dx,dy,dz} should equal centroid2 - centriod1
  Vec3 nDXYZ;
  nDXYZ.x = norm1[0]*dx;
  nDXYZ.y = norm1[1]*dy;
  nDXYZ.z = norm1[2]*dz;
  nDXYZ.abs();

  Vec3 Offset = c2 - c1;
  Offset.abs();

  Vec3 Diff = Offset - nDXYZ;
  Diff.abs();

  for (int i=0; i<3; i++) {
    if (Diff[i]>tol) return false;
  }

  // The faces are aligned across a periodic direction
  return true;
}

bool geo::comparePeriodicMPI(vector<int> &face1, vector<int> &face2)
{
  if (nDims == 2) {
    double x11, x12, y11, y12, x21, x22, y21, y22;
    x11 = xv_g(face1[0],0);  y11 = xv_g(face1[0],1);
    x12 = xv_g(face1[1],0);  y12 = xv_g(face1[1],1);
    x21 = xv_g(face2[0],0);  y21 = xv_g(face2[0],1);
    x22 = xv_g(face2[1],0);  y22 = xv_g(face2[1],1);

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

    // None of the above
    return false;
  }
  else {
    if (face1.size() != face2.size())
      return false;

    double tol = params->periodicTol;
    double dx = params->periodicDX;
    double dy = params->periodicDY;
    double dz = params->periodicDZ;

    /* --- Compare faces using normal vectors: normals should be aligned
     * and offset by norm .dot. {dx,dy,dz} --- */

    Vec3 vec1, vec2;

    // Calculate face normal & centriod for face 1
    Vec3 norm1;
    point c1;
    vec1 = point(xv_g[face1[1]]) - point(xv_g[face1[0]]);
    vec2 = point(xv_g[face1[2]]) - point(xv_g[face1[0]]);
    for (uint j=0; j<face1.size(); j++)
      c1 += xv_g[face1[j]];
    c1 /= face1.size();

    norm1[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    norm1[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    norm1[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    // Normalize
    double magNorm1 = sqrt(norm1[0]*norm1[0]+norm1[1]*norm1[1]+norm1[2]*norm1[2]);
    norm1 /= magNorm1;

    // Calculate face normal & centroid for face 2
    Vec3 norm2;
    point c2;
    vec1 = point(xv_g[face2[1]]) - point(xv_g[face2[0]]);
    vec2 = point(xv_g[face2[2]]) - point(xv_g[face2[0]]);
    for (uint j=0; j<face2.size(); j++)
      c2 += point(xv_g[face2[j]]);
    c2 /= face2.size();
    norm2[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    norm2[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    norm2[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    // Normalize
    double magNorm2 = sqrt(norm2[0]*norm2[0]+norm2[1]*norm2[1]+norm2[2]*norm2[2]);
    norm2 /= magNorm2;

    // Check for same orientation - norm1 .dot. norm2 should be +/- 1
    double dot = norm1*norm2;
    if (abs(1-abs(dot))>tol) return false;

    // Check offset distance - norm .times. {dx,dy,dz} should equal centroid2 - centriod1
    Vec3 nDXYZ;
    nDXYZ.x = norm1[0]*dx;
    nDXYZ.y = norm1[1]*dy;
    nDXYZ.z = norm1[2]*dz;
    nDXYZ.abs();

    Vec3 Offset = c2 - c1;
    Offset.abs();

    Vec3 Diff = Offset - nDXYZ;
    Diff.abs();

    for (int i=0; i<3; i++) {
      if (Diff[i]>tol) return false;
    }

    // The faces are aligned across a periodic direction
    return true;
  }
}

int geo::compareOrientation(int ic1, int f1, int ic2, int f2)
{
  if (nDims == 2) return 1;

  int nv = f2nv[c2f(ic1,f1)];

  vector<int> tmpFace1(nv), tmpFace2(nv);

  switch (ctype[ic1]) {
    case HEX:
      // Flux points arranged in 2D grid on each face oriented with each
      // dimension increasing in its +'ve direction ['btm-left' to 'top-right']
      // Node ordering reflects this: CCW from 'bottom-left' node on each face
      switch (f1) {
        case 0:
          // Bottom face  (z = -1)
          tmpFace1[0] = c2v(ic1,0);
          tmpFace1[1] = c2v(ic1,1);
          tmpFace1[2] = c2v(ic1,2);
          tmpFace1[3] = c2v(ic1,3);
          break;
        case 1:
          // Top face  (z = +1)
          tmpFace1[0] = c2v(ic1,5);
          tmpFace1[1] = c2v(ic1,4);
          tmpFace1[2] = c2v(ic1,7);
          tmpFace1[3] = c2v(ic1,6);
          break;
        case 2:
          // Left face  (x = -1)
          tmpFace1[0] = c2v(ic1,0);
          tmpFace1[1] = c2v(ic1,3);
          tmpFace1[2] = c2v(ic1,7);
          tmpFace1[3] = c2v(ic1,4);
          break;
        case 3:
          // Right face  (x = +1)
          tmpFace1[0] = c2v(ic1,2);
          tmpFace1[1] = c2v(ic1,1);
          tmpFace1[2] = c2v(ic1,5);
          tmpFace1[3] = c2v(ic1,6);
          break;
        case 4:
          // Front face  (y = -1)
          tmpFace1[0] = c2v(ic1,1);
          tmpFace1[1] = c2v(ic1,0);
          tmpFace1[2] = c2v(ic1,4);
          tmpFace1[3] = c2v(ic1,5);
          break;
        case 5:
          // Back face  (y = +1)
          tmpFace1[0] = c2v(ic1,3);
          tmpFace1[1] = c2v(ic1,2);
          tmpFace1[2] = c2v(ic1,6);
          tmpFace1[3] = c2v(ic1,7);
          break;
      }
      break;

    default:
      FatalError("Element type not supported.");
      break;
  }

  switch (ctype[ic2]) {
    case HEX:
      switch (f2) {
        case 0:
          // Bottom face  (z = -1)
          tmpFace2[0] = c2v(ic2,0);
          tmpFace2[1] = c2v(ic2,1);
          tmpFace2[2] = c2v(ic2,2);
          tmpFace2[3] = c2v(ic2,3);
          break;
        case 1:
          // Top face  (z = +1)
          tmpFace2[0] = c2v(ic2,5);
          tmpFace2[1] = c2v(ic2,4);
          tmpFace2[2] = c2v(ic2,7);
          tmpFace2[3] = c2v(ic2,6);
          break;
        case 2:
          // Left face  (x = -1)
          tmpFace2[0] = c2v(ic2,0);
          tmpFace2[1] = c2v(ic2,3);
          tmpFace2[2] = c2v(ic2,7);
          tmpFace2[3] = c2v(ic2,4);
          break;
        case 3:
          // Right face  (x = +1)
          tmpFace2[0] = c2v(ic2,2);
          tmpFace2[1] = c2v(ic2,1);
          tmpFace2[2] = c2v(ic2,5);
          tmpFace2[3] = c2v(ic2,6);
          break;
        case 4:
          // Front face  (y = -1)
          tmpFace2[0] = c2v(ic2,1);
          tmpFace2[1] = c2v(ic2,0);
          tmpFace2[2] = c2v(ic2,4);
          tmpFace2[3] = c2v(ic2,5);
          break;
        case 5:
          // Back face  (y = +1)
          tmpFace2[0] = c2v(ic2,3);
          tmpFace2[1] = c2v(ic2,2);
          tmpFace2[2] = c2v(ic2,6);
          tmpFace2[3] = c2v(ic2,7);
          break;
      }
      break;

    default:
      FatalError("Element type not supported.");
      break;
  }

  // Now, compare the two faces to see the relative orientation [rotation]
  if      (tmpFace1[0] == tmpFace2[0]) return 0;
  else if (tmpFace1[1] == tmpFace2[0]) return 1;
  else if (tmpFace1[2] == tmpFace2[0]) return 2;
  else if (tmpFace1[3] == tmpFace2[0]) return 3;
  else if (checkPeriodicFaces3D(tmpFace1,tmpFace2)) {
    point c1, c2;
    for (auto iv:tmpFace1) c1 += point(xv[iv]);
    for (auto iv:tmpFace2) c2 += point(xv[iv]);
    c1 /= tmpFace1.size();
    c2 /= tmpFace2.size();
    Vec3 fDist = c2 - c1;
    fDist /= sqrt(fDist*fDist); // Normalize

    point pt1;
    point pt2 = point(xv[tmpFace2[0]]);

    for (int i=0; i<4; i++) {
      pt1 = point(xv[tmpFace1[i]]);
      Vec3 ptDist = pt2 - pt1;        // Vector between points
      ptDist /= sqrt(ptDist*ptDist); // Normalize

      double dot = fDist*ptDist;
      if (abs(1-abs(dot))<params->periodicTol) return i; // These points align
    }
    // Matching points not found
    FatalError("Unable to orient periodic faces using simple algorithm.");
  }
  else FatalError("Internal faces improperly matched.");

}

int geo::compareOrientationMPI(int ic1, int f1, int ic2, int f2, int isPeriodic=0)
{
  if (nDims == 2) return 1;

  int nv = f2nv[c2f(ic1,f1)];

  vector<int> tmpFace1(nv), tmpFace2(nv);

  int icg = ic2icg[ic1];
  switch (ctype[ic1]) {
    case HEX:
      // Flux points arranged in 2D grid on each face oriented with each
      // dimension increasing in its +'ve direction ['btm-left' to 'top-right']
      // Node ordering reflects this: CCW from 'bottom-left' node on each face
      switch (f1) {
        case 0:
          // Bottom face  (z = -1)
          tmpFace1[0] = c2v_g(icg,0);
          tmpFace1[1] = c2v_g(icg,1);
          tmpFace1[2] = c2v_g(icg,2);
          tmpFace1[3] = c2v_g(icg,3);
          break;
        case 1:
          // Top face  (z = +1)
          tmpFace1[0] = c2v_g(icg,5);
          tmpFace1[1] = c2v_g(icg,4);
          tmpFace1[2] = c2v_g(icg,7);
          tmpFace1[3] = c2v_g(icg,6);
          break;
        case 2:
          // Left face  (x = -1)
          tmpFace1[0] = c2v_g(icg,0);
          tmpFace1[1] = c2v_g(icg,3);
          tmpFace1[2] = c2v_g(icg,7);
          tmpFace1[3] = c2v_g(icg,4);
          break;
        case 3:
          // Right face  (x = +1)
          tmpFace1[0] = c2v_g(icg,2);
          tmpFace1[1] = c2v_g(icg,1);
          tmpFace1[2] = c2v_g(icg,5);
          tmpFace1[3] = c2v_g(icg,6);
          break;
        case 4:
          // Front face  (y = -1)
          tmpFace1[0] = c2v_g(icg,1);
          tmpFace1[1] = c2v_g(icg,0);
          tmpFace1[2] = c2v_g(icg,4);
          tmpFace1[3] = c2v_g(icg,5);
          break;
        case 5:
          // Back face  (y = +1)
          tmpFace1[0] = c2v_g(icg,3);
          tmpFace1[1] = c2v_g(icg,2);
          tmpFace1[2] = c2v_g(icg,6);
          tmpFace1[3] = c2v_g(icg,7);
          break;
      }
      break;

    default:
      FatalError("Element type not supported.");
      break;
  }

  switch (ctype_g[ic2]) {
    case HEX:
      switch (f2) {
        case 0:
          // Bottom face  (z = -1)
          tmpFace2[0] = c2v_g(ic2,0);
          tmpFace2[1] = c2v_g(ic2,1);
          tmpFace2[2] = c2v_g(ic2,2);
          tmpFace2[3] = c2v_g(ic2,3);
          break;
        case 1:
          // Top face  (z = +1)
          tmpFace2[0] = c2v_g(ic2,5);
          tmpFace2[1] = c2v_g(ic2,4);
          tmpFace2[2] = c2v_g(ic2,7);
          tmpFace2[3] = c2v_g(ic2,6);
          break;
        case 2:
          // Left face  (x = -1)
          tmpFace2[0] = c2v_g(ic2,0);
          tmpFace2[1] = c2v_g(ic2,3);
          tmpFace2[2] = c2v_g(ic2,7);
          tmpFace2[3] = c2v_g(ic2,4);
          break;
        case 3:
          // Right face  (x = +1)
          tmpFace2[0] = c2v_g(ic2,2);
          tmpFace2[1] = c2v_g(ic2,1);
          tmpFace2[2] = c2v_g(ic2,5);
          tmpFace2[3] = c2v_g(ic2,6);
          break;
        case 4:
          // Front face  (y = -1)
          tmpFace2[0] = c2v_g(ic2,1);
          tmpFace2[1] = c2v_g(ic2,0);
          tmpFace2[2] = c2v_g(ic2,4);
          tmpFace2[3] = c2v_g(ic2,5);
          break;
        case 5:
          // Back face  (y = +1)
          tmpFace2[0] = c2v_g(ic2,3);
          tmpFace2[1] = c2v_g(ic2,2);
          tmpFace2[2] = c2v_g(ic2,6);
          tmpFace2[3] = c2v_g(ic2,7);
          break;
      }
      break;

    default:
      FatalError("Element type not supported.");
      break;
  }

  // Now, compare the two faces to see the relative orientation [rotation]
  if      (tmpFace1[0] == tmpFace2[0]) return 0;
  else if (tmpFace1[1] == tmpFace2[0]) return 1;
  else if (tmpFace1[2] == tmpFace2[0]) return 2;
  else if (tmpFace1[3] == tmpFace2[0]) return 3;
  else if (isPeriodic) {
    if (!comparePeriodicMPI(tmpFace1,tmpFace2))
      FatalError("Periodic MPI faces improperly matched.");

    point c1, c2;
    for (auto iv:tmpFace1) c1 += point(xv_g[iv]);
    for (auto iv:tmpFace2) c2 += point(xv_g[iv]);
    c1 /= tmpFace1.size();
    c2 /= tmpFace2.size();
    Vec3 fDist = c2 - c1;
    fDist /= sqrt(fDist*fDist); // Normalize

    point pt1;
    point pt2 = point(xv_g[tmpFace2[0]]);

    for (int i=0; i<4; i++) {
      pt1 = point(xv_g[tmpFace1[i]]);
      Vec3 ptDist = pt2 - pt1;        // Vector between points
      ptDist /= sqrt(ptDist*ptDist); // Normalize

      double dot = fDist*ptDist;
      if (abs(1-abs(dot))<params->periodicTol) return i; // These points align
    }
    // Matching points not found
    FatalError("Unable to orient periodic MPI faces using simple algorithm.");
  }
  else FatalError("MPI faces improperly matched.");

}

void geo::partitionMesh(void)
{
#ifndef _NO_MPI

  if (nproc <= 1) {
    gridComm = MPI_COMM_WORLD;
    return;
  }

  if (meshType == OVERSET_MESH) {
    // Partitioning each grid independantly; local 'grid rank' is the important rank
    rank = gridRank;
    nproc = nProcGrid;

    if (nproc <= 1) return; // No additional partitioning needed

    if (rank == 0) cout << "Geo: Partitioning mesh block " << gridID << " across " << nProcGrid << " processes" << endl;
    if (rank == 0) cout << "Geo:   Number of elements in block " << gridID << " : " << nEles << endl;
  }
  else {
    if (rank == 0) cout << "Geo: Partitioning mesh across " << nproc << " processes" << endl;
    if (rank == 0) cout << "Geo:   Number of elements globally: " << nEles << endl;

    gridComm = MPI_COMM_WORLD;
  }

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
  epart.resize(nEles);
  vector<int> npart(nVerts);

  // int errVal = METIS PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize,
  // idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval,idx_t *epart, idx_t *npart)

  int ncommon; // 2 for 2D, ~3 for 3D [#nodes per face: 2 for quad/tri, 3 for tet, 4 for hex]
  if (nDims == 2) ncommon = 2;
  else if (nDims == 3) ncommon = 4;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // Reduces # of final MPI faces
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
  options[METIS_OPTION_NCUTS] = 5;  // Allows better partitioning (less cuts) to be found [at negligible expense for CFD grids]

  METIS_PartMeshDual(&nEles,&nVerts,eptr.data(),eind.data(),NULL,NULL,
                     &ncommon,&nproc,NULL,options,&objval,epart.data(),npart.data());

  // Copy data to the global arrays & reset local arrays
  nEles_g   = nEles;
  nVerts_g  = nVerts;
  c2v_g     = c2v;      c2v.setup(0,0);
  xv_g      = xv;       xv.setup(0,0);
  ctype_g   = ctype;    ctype.resize(0);
  c2nv_g    = c2nv;     c2nv.resize(0);
  c2ne_g    = c2nf;     c2nf.resize(0);
  bndPts_g  = bndPts;   bndPts.setup(0,0);
  nBndPts_g = nBndPts;  nBndPts.resize(nBounds);

  // Each processor will now grab its own data according to its rank (proc ID)
  for (int i=0; i<nEles; i++) {
    if (epart[i] == rank) {
      c2v.insertRow(c2v_g[i],-1,c2nv_g[i]);
      ic2icg.push_back(i);
      ctype.push_back(ctype_g[i]);
      c2nv.push_back(c2nv_g[i]);
      c2nf.push_back(c2ne_g[i]);
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
  vector<int> ivg2iv;
  ivg2iv.assign(nVerts_g,-1);

  // Transfer over all needed vertices to local array
  int nv = 0;
  for (auto &iv: myNodes) {
    xv.insertRow(xv_g.getRow(iv));
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

  // Lastly, update c2v from global --> local node IDs
  std::transform(c2v.getData(),c2v.getData()+c2v.getSize(),c2v.getData(), [=](int ivg){return ivg2iv[ivg];});

  if (meshType == OVERSET_MESH)
    cout << "Geo:   Grid " << gridID << " on rank " << rank << ": nEles = " << nEles << endl;
  else
    cout << "Geo:   On rank " << rank << ": nEles = " << nEles << endl;

  if (rank == 0) cout << "Geo: Done partitioning mesh" << endl;
#endif
}

void geo::getMpiPartitions(void)
{
#ifndef _NO_MPI

  if (nproc <= 1) {
    gridComm = MPI_COMM_WORLD;
    return;
  }

  if (meshType == OVERSET_MESH) {
    // Partitioning each grid independantly; local 'grid rank' is the important rank
    rank = gridRank;
    nproc = nProcGrid;

    if (rank == 0) cout << "Geo: Partitioning mesh block " << gridID << " across " << nProcGrid << " processes" << endl;
    if (rank == 0) cout << "Geo:   Number of elements in block " << gridID << " : " << nEles << endl;

    if (nproc <= 1) return; // No additional partitioning needed
  }
  else {
    if (rank == 0) cout << "Geo: Partitioning mesh across " << nproc << " processes" << endl;
    if (rank == 0) cout << "Geo:   Number of elements globally: " << nEles << endl;

    gridComm = MPI_COMM_WORLD;
  }

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
  epart.resize(nEles);
  vector<int> npart(nVerts);

  // int errVal = METIS PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize,
  // idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval,idx_t *epart, idx_t *npart)

  int ncommon; // 2 for 2D, ~3 for 3D [#nodes per face: 2 for quad/tri, 3 for tet, 4 for hex]
  if (nDims == 2) ncommon = 2;
  else if (nDims == 3) ncommon = 4;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // Reduces # of final MPI faces
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
  options[METIS_OPTION_NCUTS] = 5;  // Allows better partitioning (less cuts) to be found [at negligible expense for CFD grids]

  METIS_PartMeshDual(&nEles,&nVerts,eptr.data(),eind.data(),NULL,NULL,
                     &ncommon,&nproc,NULL,options,&objval,epart.data(),npart.data());
#endif
}

void geo::partitionFromEpart(const vector<int>& _epart)
{
#ifndef _NO_MPI
  epart = _epart;
  if (params->nproc <= 1) {
    gridComm = MPI_COMM_WORLD;
    return;
  }

  if (meshType == OVERSET_MESH) {
    // Partitioning each grid independantly; local 'grid rank' is the important rank
    rank = gridRank;
    nproc = nProcGrid;
    if (nproc <= 1) return; // No additional partitioning needed
  }
  else {
    gridComm = MPI_COMM_WORLD;
  }

  // Copy data to the global arrays & reset local arrays
  nEles_g   = nEles;
  nVerts_g  = nVerts;
  c2v_g     = c2v;      c2v.setup(0,0);
  xv_g      = xv;       xv.setup(0,0);
  ctype_g   = ctype;    ctype.resize(0);
  c2nv_g    = c2nv;     c2nv.resize(0);
  c2ne_g    = c2nf;     c2nf.resize(0);
  bndPts_g  = bndPts;   bndPts.setup(0,0);
  nBndPts_g = nBndPts;  nBndPts.resize(nBounds);

  // Each processor will now grab its own data according to its rank (proc ID)
  for (int i=0; i<nEles; i++) {
    if (epart[i] == rank) {
      c2v.insertRow(c2v_g[i],-1,c2nv_g[i]);
      ic2icg.push_back(i);
      ctype.push_back(ctype_g[i]);
      c2nv.push_back(c2nv_g[i]);
      c2nf.push_back(c2ne_g[i]);
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
  vector<int> ivg2iv;
  ivg2iv.assign(nVerts_g,-1);

  // Transfer over all needed vertices to local array
  int nv = 0;
  for (auto &iv: myNodes) {
    xv.insertRow(xv_g.getRow(iv));
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

  // Lastly, update c2v from global --> local node IDs
  std::transform(c2v.getData(),c2v.getData()+c2v.getSize(),c2v.getData(), [=](int ivg){return ivg2iv[ivg];});

  if (meshType == OVERSET_MESH)
    cout << "Geo:   Grid " << gridID << " on rank " << rank << ": nEles = " << nEles << endl;
  else
    cout << "Geo:   On rank " << rank << ": nEles = " << nEles << endl;

  if (rank == 0) cout << "Geo: Done partitioning mesh" << endl;
#endif
}

void geo::moveMesh(double rkVal)
{
  double rkTime = params->time + params->dt*rkVal;
  params->rkTime = rkTime;

  switch (params->motion) {
    case 1: {
      #pragma omp parallel for
      for (int iv=0; iv<nVerts; iv++) {
        /// Taken from Kui, AIAA-2010-5031-661
        xv(iv,0) = xv0[iv].x + 2*sin(pi*xv0[iv].x/10.)*sin(pi*xv0[iv].y/10.)*sin(2*pi*rkTime/10.);
        xv(iv,1) = xv0[iv].y + 2*sin(pi*xv0[iv].x/10.)*sin(pi*xv0[iv].y/10.)*sin(2*pi*rkTime/10.);
        gridVel(iv,0) = 4.*pi/10.*sin(pi*xv0[iv].x/10.)*sin(pi*xv0[iv].y/10.)*cos(2*pi*rkTime/10.);
        gridVel(iv,1) = 4.*pi/10.*sin(pi*xv0[iv].x/10.)*sin(pi*xv0[iv].y/10.)*cos(2*pi*rkTime/10.);
      }
      break;
    }
    case 2: {
      double t0 = 10.*sqrt(5.);
      if (nDims == 2) {
        #pragma omp parallel for
        for (int iv=0; iv<nVerts; iv++) {
          /// Taken from Liang-Miyaji
          xv(iv,0) = xv0[iv].x + sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(4*pi*rkTime/t0);
          xv(iv,1) = xv0[iv].y + sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(8*pi*rkTime/t0);
          gridVel(iv,0) = 4.*pi/t0*sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*cos(4*pi*rkTime/t0);
          gridVel(iv,1) = 8.*pi/t0*sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*cos(8*pi*rkTime/t0);
        }
      }
      else {
        #pragma omp parallel for
        for (int iv=0; iv<nVerts; iv++) {
          /// Taken from Liang-Miyaji
          xv(iv,0) = xv0[iv].x + sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*sin(4*pi*rkTime/t0);
          xv(iv,1) = xv0[iv].y + sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*sin(8*pi*rkTime/t0);
          xv(iv,2) = xv0[iv].z + sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*sin(4*pi*rkTime/t0);
          gridVel(iv,0) = 4.*pi/t0*sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*cos(4*pi*rkTime/t0);
          gridVel(iv,1) = 8.*pi/t0*sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*cos(8*pi*rkTime/t0);
          gridVel(iv,2) = 4.*pi/t0*sin(pi*xv0[iv].x/5.)*sin(pi*xv0[iv].y/5.)*sin(pi*xv0[iv].z/5.)*cos(4*pi*rkTime/t0);
        }
      }
      break;
    }
    case 3: {
      if (gridID==0) {
        /// Liangi-Miyaji with easily-modifiable domain width
        double t0 = 10.*sqrt(5.);
        double width = 5.;
        #pragma omp parallel for
        for (int iv=0; iv<nVerts; iv++) {
          xv(iv,0) = xv0[iv].x + sin(pi*xv0[iv].x/width)*sin(pi*xv0[iv].y/width)*sin(4*pi*rkTime/t0);
          xv(iv,1) = xv0[iv].y + sin(pi*xv0[iv].x/width)*sin(pi*xv0[iv].y/width)*sin(8*pi*rkTime/t0);
          gridVel(iv,0) = 4.*pi/t0*sin(pi*xv0[iv].x/width)*sin(pi*xv0[iv].y/width)*cos(4*pi*rkTime/t0);
          gridVel(iv,1) = 8.*pi/t0*sin(pi*xv0[iv].x/width)*sin(pi*xv0[iv].y/width)*cos(8*pi*rkTime/t0);
        }
      }
      break;
    }
    case 4: {
      /// Rigid oscillation in a circle
      if (gridID==0) {
        double Ax = params->moveAx; // Amplitude  (m)
        double Ay = params->moveAy; // Amplitude  (m)
        double fx = params->moveFx; // Frequency  (Hz)
        double fy = params->moveFy; // Frequency  (Hz)
        for (int iv=0; iv<nVerts; iv++) {
          xv(iv,0) = xv0[iv].x + Ax*sin(2.*pi*fx*rkTime);
          //xv(iv,1) = xv0[iv].y + Ay*sin(2.*pi*fy*rkTime);
          xv(iv,1) = xv0[iv].y + Ay*(1-cos(2.*pi*fy*rkTime));
          gridVel(iv,0) = 2.*pi*fx*Ax*cos(2.*pi*fx*rkTime);
          //gridVel(iv,1) = 2.*pi*fy*Ay*cos(2.*pi*fy*rkTime);
          gridVel(iv,1) = 2.*pi*fy*Ay*sin(2.*pi*fy*rkTime);
        }
      }
    }
  }
}
