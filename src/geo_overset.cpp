/*!
 * \file geo_overset.cpp
 * \brief Overset-related methods for the geo class
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */

#include "geo.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>

#include "ele.hpp"
#include "face.hpp"
#include "overFace.hpp"
#include "overComm.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#include "metis.h"
#endif

void geo::splitGridProcs(void)
{
#ifndef _NO_MPI
  // Split the processes among the overset grids such that they are roughly balanced

  /* --- Read Number of Elements in Each Grid --- */

  vector<int> nElesGrid(nGrids);
  int nElesTotal = 0;

  for (int i=0; i<nGrids; i++) {
    ifstream meshFile;
    string str;
    string fileName = params->oversetGrids[i];

    meshFile.open(fileName.c_str());
    if (!meshFile.is_open())
      FatalError("Unable to open mesh file.");

    // Move cursor to $Elements
    meshFile.clear();
    meshFile.seekg(0, ios::beg);
    while(1) {
      getline(meshFile,str);
      if (str.find("$Elements")!=string::npos) break;
      if(meshFile.eof()) FatalError("$Elements tag not found in Gmsh file!");
    }

    // Read total number of interior + boundary elements
    meshFile >> nElesGrid[i];
    meshFile.close();

    nElesTotal += nElesGrid[i];
  }

  /* --- Balance the processes across the grids --- */

  nProcsGrid.resize(nGrids);
  for (int i=0; i<nGrids; i++) {
    double eleRatio = (double)nElesGrid[i]/nElesTotal;
    nProcsGrid[i] = round(eleRatio*nproc);
  }

  /* --- Get the final gridID for this rank --- */

  int g = 0;
  int procSum = nProcsGrid[0];
  while (procSum < rank+1 && g<nGrids-1) {
    g++;
    procSum += nProcsGrid[g];
  }
  gridID = g;

  /* --- Split MPI Processes Based Upon gridID: Create MPI_Comm for each grid --- */

  MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);

  MPI_Comm_rank(gridComm,&gridRank);
  MPI_Comm_size(gridComm,&nProcGrid);

  cout << "rank " << rank << ", gridRank " << gridRank << ", gridID " << gridID << endl;

  gridIdList.resize(nproc);
  MPI_Allgather(&gridID,1,MPI_INT,gridIdList.data(),1,MPI_INT,MPI_COMM_WORLD);
#endif
}

void geo::registerGridDataTIOGA(void)
{
#ifndef _NO_MPI
  /* Note that this function should only be needed once during preprocessing */

  if (gridRank == 0)
    cout << "Geo: Grid " << gridID << ": Setting up TIOGA interface" << endl;

  // Allocate TIOGA grid processor
  tg = make_shared<tioga>();

  tg->setCommunicator(MPI_COMM_WORLD,params->rank,params->nproc);

  // Setup iwall, iover (nodes on wall & overset boundaries)
  iover.resize(0);
  iwall.resize(0);
  nodeType.assign(nVerts,NORMAL_NODE);
  for (int ib=0; ib<nBounds; ib++) {
    if (bcList[ib] == OVERSET) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iover.push_back(bndPts(ib,iv));
        nodeType[bndPts(ib,iv)] = OVERSET_NODE;
      }
    }
    else if (bcList[ib] == SLIP_WALL || bcList[ib] == ADIABATIC_NOSLIP || bcList[ib] == ISOTHERMAL_NOSLIP) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iwall.push_back(bndPts(ib,iv));
        nodeType[bndPts(ib,iv)] = BOUNDARY_NODE;
      }
    }
    else {
      nodeType[bndPts(ib,ib)] = BOUNDARY_NODE;
    }
  }

  int nwall = iwall.size();
  int nover = iover.size();
  int ntypes = 1;           //! Number of element types in grid block
  nodesPerCell = new int[1];
  nodesPerCell[0] = 8;      //! Number of nodes per element for each element type (but only one type so far)
  iblank.resize(nVerts);
  iblankCell.resize(nEles);
  iblankFace.resize(nFaces);

  if (c2v.getDim1() > 8) {
    // Quadratic elements present; setup a different c2v specifically for TIOGA
    tg_c2v.setup(nEles,8);
    for (int ic=0; ic<nEles; ic++) {
      for (int j=0; j<8; j++) {
        tg_c2v(ic,j) = c2v(ic,j);
      }
    }
    // Need an int**, even if only have one element type
    conn[0] = tg_c2v.getData();
  }
  else {
    // Need an int**, even if only have one element type
    conn[0] = c2v.getData();
  }

  tg->registerGridData(gridID,nVerts,xv.getData(),iblank.data(),nwall,nover,iwall.data(),
                       iover.data(),ntypes,nodesPerCell,&nEles,&conn[0]);

  // Get a list of all cells which have an overset-boundary face (for later use with blanking)
  for (int i=0; i<nBndFaces; i++) {
    if (bcType[i]==OVERSET) {
      overCells.insert(f2c(bndFaces[i],0));
    }
  }
#endif
}

void geo::setupOverset2D(void)
{
#ifndef _NO_MPI
  OComm = make_shared<overComm>();

  OComm->setup(params,nGrids,gridID,gridRank,nProcGrid,gridIdList);

  setNodeTypes2D();

  iblank.resize(nVerts);
  iblankCell.resize(nEles);
  iblankFace.resize(nFaces);

  OComm->setIblanks2D(xv,overFaceNodes,wallFaceNodes,iblank);

  // Now use new nodal iblanks to set cell and face iblanks
  setCellIblanks();
#endif
}

void geo::setNodeTypes2D(void)
{
  // Set node types for 'mandatory' blanking values, and get wall & overset nodes
  iover.resize(0);
  iwall.resize(0);
  nodeType.assign(nVerts,NORMAL_NODE);
  for (int ib=0; ib<nBounds; ib++) {
    if (bcList[ib] == OVERSET) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iover.push_back(bndPts(ib,iv));
        nodeType[bndPts(ib,iv)] = OVERSET_NODE;
      }
    }
    else if (bcList[ib] == SLIP_WALL || bcList[ib] == ADIABATIC_NOSLIP || bcList[ib] == ISOTHERMAL_NOSLIP) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iwall.push_back(bndPts(ib,iv));
        nodeType[bndPts(ib,iv)] = BOUNDARY_NODE;
      }
    }
    else {
      nodeType[bndPts(ib,ib)] = BOUNDARY_NODE;
    }
  }

  overFaceNodes.setup(0,0);
  wallFaceNodes.setup(0,0);
  for (int bf=0; bf<nBndFaces; bf++) {
    if (bcType[bf] == SLIP_WALL || bcType[bf] == ADIABATIC_NOSLIP || bcType[bf] == ISOTHERMAL_NOSLIP) {
      int ff = bndFaces[bf];
      wallFaceNodes.insertRow({f2v(ff,0),f2v(ff,1)});
    }
    else if (bcType[bf] == OVERSET) {
      int ff = bndFaces[bf];
      overFaceNodes.insertRow({f2v(ff,0),f2v(ff,1)});
    }
  }
}

void geo::updateADT(void)
{
#ifndef _NO_MPI
  if (nDims == 3) {
    // Pre-process the grids
    tg->profile();

    // Have TIOGA perform the nodal overset connectivity (set nodal iblanks)
    tg->performConnectivity();
  }
#endif
}

void geo::updateBlanking(void)
{
#ifndef _NO_MPI
  updateADT();

  if (nDims == 2) {
    // TIOGA only for 3D, so use my own 2D hole-cutting implementation
    OComm->setIblanks2D(xv,overFaceNodes,wallFaceNodes,iblank);
  }

  // Now use new nodal iblanks to set cell and face iblanks
  setCellIblanks();
#endif
}

void geo::setCellIblanks(void)
{
  // Use the TIOGA-supplied nodal iblank values, set iblank values for all cells and faces

  // Only needed for moving grids: List of current hole cells
  holeCells.clear();
  blankCells.clear();
  unblankCells.clear();
  for (int ic=0; ic<nEles; ic++)
    if (iblankCell[ic] == HOLE)
      holeCells.insert(ic);

  iblankCell.assign(nEles,NORMAL);

  // First, blank any fringe vertices which should be treated as hole vertices
  vector<int> iblank1(nVerts,HOLE);
  for (int iv=0; iv<nVerts; iv++) {
    if (iblank[iv] == FRINGE) {
      for (int j=0; j<v2nv[iv]; j++) {
        if (iblank[v2v(iv,j)] == NORMAL) {
          iblank1[iv] = NORMAL;
        }
      }
    }
  }

  for (int iv=0; iv<nVerts; iv++) {
    if (iblank1[iv] == NORMAL) iblank[iv] = NORMAL;
  }

  for (int iv=0; iv<nVerts; iv++) {
    if (iblank[iv] == FRINGE) {
      int nfringe = 0;
      for (int j=0; j<v2nv[iv]; j++) {
        if ((iblank[v2v(iv,j)] == FRINGE && nodeType[v2v(iv,j)] == NORMAL_NODE) || iblank[v2v(iv,j)] == HOLE ) {
          nfringe++;
        }
      }
      if (nfringe == v2nv[iv])
        iblank[iv] = HOLE; // HOLE
    }
  }

  // Next, blank all cells which contain a hole node

  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nv[ic]; j++) {
      int iv = c2v(ic,j);
      if (iblank[iv] == HOLE) {
        iblankCell[ic] = HOLE;

        // Only needed for moving grids: Existing cells which must be removed from solver
        if (!holeCells.count(ic))
          blankCells.insert(ic);

        break;
      }
    }
  }

  // Only needed for moving grids: Get cells which  must be 'un-blanked'
  for (auto &ic:holeCells)
    if (iblankCell[ic] == NORMAL)
      unblankCells.insert(ic);
}

void geo::setFaceIblanks(void)
{
  iblankFace.assign(nFaces,NORMAL);

  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic]!=HOLE) continue;

    for (int j=0; j<c2nf[ic]; j++) {
      if (c2c(ic,j)>=0) {
        // Internal face
        if (iblankCell[c2c(ic,j)] == HOLE) {
          iblankFace[c2f(ic,j)] = HOLE;
        } else if (iblankCell[c2c(ic,j)] == NORMAL) {
          iblankFace[c2f(ic,j)] = FRINGE;
        }
      } else {
        // Boundary or MPI face
        iblankFace[c2f(ic,j)] = HOLE;
      }
    }
  }
}

void geo::matchOversetDonors(vector<shared_ptr<ele>> &eles, vector<superMesh> &donors)
{
#ifndef _NO_MPI

#endif
}

void geo::processBlanks(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces)
{
#ifndef _NO_MPI
  /* --- Set blank/unblank faces for all elements to be blanked --- */

  set<int> blankIFaces, blankMFaces, blankOFaces, ubOFaces;
  for (auto &ic:blankCells) {
    for (int j=0; j<c2nf[ic]; j++) {
      int ftype = currFaceType[c2f(ic,j)];
      switch (ftype) {
        case INTERNAL:
          blankIFaces.insert(c2f(ic,j));
          break;
        case BOUNDARY:
          blankIFaces.insert(c2f(ic,j));
          break;
        case MPI_FACE:
          blankMFaces.insert(c2f(ic,j));
          break;
        case OVER_FACE:
          blankOFaces.insert(c2f(ic,j));
          break;
        case HOLE_FACE:
          //FatalError("Face of blankCell is already blanked.");
          break;
        default:
          FatalError("Face of blankCell has unknown currFaceType.");
          break;
      }
    }
  }

  // Find all the internal faces which must be replaced with overset faces
  for (auto &ff:blankIFaces) {
    int ic1 = f2c(ff,0);
    int ic2 = f2c(ff,1);
    if ( (iblankCell[ic1]==NORMAL && iblankCell[ic2]==HOLE) ||
         (iblankCell[ic1]==HOLE && iblankCell[ic2]==NORMAL) )
      ubOFaces.insert(ff);
  }

  // Figure out whether any other MPI faces must be replaced with overset faces
  vector<int> blankMpi;
  for (auto &ff:blankMFaces) blankMpi.push_back(ff);

  int nblankMpi = blankMpi.size();
  vector<int> nblankMpi_rank(nProcGrid);
  MPI_Allgather(&nblankMpi,1,MPI_INT,nblankMpi_rank.data(),1,MPI_INT,gridComm);

  int sum = 0;
  vector<int> recvCnts(nProcGrid);
  vector<int> recvDisp(nProcGrid);
  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nblankMpi_rank[i];
    if (i>0)
      recvDisp[i] = recvDisp[i-1]+recvCnts[i-1];
    sum += recvCnts[i];
  }
  vector<int> blankMpi_rank(sum);
  MPI_Allgatherv(blankMpi.data(),blankMpi.size(),MPI_INT,blankMpi_rank.data(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

  for (int F=0; F<nMpiFaces; F++) {
    int ff = mpiFaces[F];
    int p = procR[F];
    int f2 = faceID_R[F];

    bool isBlanked = false;
    for (int j=0; j<recvCnts[p]; j++) {
      if (blankMpi_rank[recvDisp[p]+j]==f2) {
        isBlanked = true;
        break;
      }
    }

    // If MPI face is blanked on other side but not this one, must replace with an overFace
    if (isBlanked && !blankMFaces.count(ff)) {
      blankMFaces.insert(ff);
      ubOFaces.insert(ff);
    }
  }

  set<int> ubIFaces, ubMFaces;

  removeEles(eles,blankCells);
  removeFaces(faces,mFaces,oFaces,blankIFaces,blankMFaces,blankOFaces);
  insertFaces(eles,faces,mFaces,oFaces,ubIFaces,ubMFaces,ubOFaces);

  for (auto &iface:faces) {
    iface->getPointers();
    iface->getPointersRight();
  }
  for (auto &mface:mFaces) mface->getPointers();
  for (auto &oface:oFaces) oface->getPointers();
#endif
}

void geo::processUnblanks(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces)
{
  /* --- Set Unblank/Blank Faces for All Unblank Elements --- */

  set<int> ubIntFaces, ubMpiFaces, ubOFaces;
  set<int> blankIFaces, blankMFaces, blankOFaces;
  for (auto &ic:unblankCells) {
    for (int j=0; j<c2nf[ic]; j++) {
      int ic2 = c2c(ic,j);
      int ff2 = c2f(ic,j);
      if (ic2>=0) {
        if (iblankCell[ic2]==NORMAL || unblankCells.count(ic2))
          ubIntFaces.insert(ff2);
        else
          ubOFaces.insert(ff2);
      } else {
        // Boundary or MPI face
        if (faceType[ff2]==MPI_FACE)
          ubMpiFaces.insert(ff2);
        else
          ubIntFaces.insert(ff2);
      }
    }
  }

#ifndef _NO_MPI
  // Figure out whether MPI faces are to be unblanked as MPI or as overset
  vector<int> ubMpi;
  for (auto &ff:ubMpiFaces) ubMpi.push_back(ff);

  int nUbMPI = ubMpi.size();
  vector<int> nubMpi_rank(nProcGrid);
  MPI_Allgather(&nUbMPI,1,MPI_INT,nubMpi_rank.data(),1,MPI_INT,gridComm);

  int sum = 0;
  vector<int> recvCnts(nProcGrid);
  vector<int> recvDisp(nProcGrid);
  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nubMpi_rank[i];
    if (i>0)
      recvDisp[i] = recvDisp[i-1]+recvCnts[i-1];
    sum += recvCnts[i];
  }
  vector<int> ubMpi_rank(sum);
  MPI_Allgatherv(ubMpi.data(),ubMpi.size(),MPI_INT,ubMpi_rank.data(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

  for (int F=0; F<nMpiFaces; F++) {
    int ff = mpiFaces[F];
    if (ubMpiFaces.count(ff) || overFaces.count(ff)) {
      int p = procR[F];
      int f2 = faceID_R[F];

      bool isUnblanked = false;
      for (int j=0; j<recvCnts[p]; j++) {
        if (ubMpi_rank[recvDisp[p]+j]==f2) {
          isUnblanked = true;
          break;
        }
      }

      if (ubMpiFaces.count(ff) && !isUnblanked) {
        // If MPI face is not unblanked on other rank, move to overFaces
        ubMpiFaces.erase(ff);
        ubOFaces.insert(ff);
      } else if (currFaceType[ff]==OVER_FACE && isUnblanked) {
        // Current overset face must be replaced by MPI face
        blankOFaces.insert(ff);
        ubMpiFaces.insert(ff);
      }
    }
  }
#endif

  // Now, figure out what faces must be removed due to being replaced by other type
  // For cell unblanking, the only possibility for face blanking is overset faces

  for (auto &ff:ubIntFaces)
    if (overFaces.count(ff))
      blankOFaces.insert(ff);

  for (auto &ff:ubMpiFaces)
    if (overFaces.count(ff))
      blankOFaces.insert(ff);

  insertEles(eles,unblankCells);
  removeFaces(faces,mFaces,oFaces,blankIFaces,blankMFaces,blankOFaces);
  insertFaces(eles,faces,mFaces,oFaces,ubIntFaces,ubMpiFaces,ubOFaces);

  for (auto &iface:faces) {
    iface->getPointers();
    iface->getPointersRight();
  }
  for (auto &mface:mFaces) mface->getPointers();
  for (auto &oface:oFaces) oface->getPointers();
}

void geo::removeEles(vector<shared_ptr<ele>> &eles, set<int> &blankEles)
{
  /* --- Remove Newly-Blanked Elements --- */

  for (auto &ic:blankEles) {
    if (ic<0) continue;
    int ind = eleMap[ic];
    if (ind<0) FatalError("Should not have marked a hole cell for blanking!");
    eles.erase(eles.begin()+ind,eles.begin()+ind+1);
    eleMap[ic] = -1;

    // Update the map
    for (int k=ic+1; k<nEles; k++)
      if (eleMap[k]>=0)
        eleMap[k]--;
  }
}

void geo::insertEles(vector<shared_ptr<ele>> &eles, set<int> &ubEles)
{
  /* --- Setup & Insert Unblanked Elements --- */

  for (auto &ic:ubEles) {
    // Find the next-lowest index
    int ind = eleMap[ic];
    if (ind>=0) FatalError("Should not have marked a non-hole cell for un-blanking! Is eleMap wrong?");
    int j = 1;
    while (ind < 0 && j<=ic) {
      ind = eleMap[ic-j];
      j++;
    }
    ind++;

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

    e->setup(params,this);

    eles.insert(eles.begin()+ind,1,e);

    // Update the map
    eleMap[ic] = ind;
    for (int k=ic+1; k<nEles; k++)
      if (eleMap[k]>=0)
        eleMap[k]++;
  }
}

void geo::insertFaces(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces,
                      set<int> &ubIFaces, set<int> &ubMFaces, set<int> &ubOFaces)
{
  /* --- Create Newly-Unblanked Faces Corresponding to Unblanked Eles --- */

  for (auto &ff:ubIFaces) {
    if (ff<0) continue;

    if (faceType[ff] == INTERNAL)
    {
      shared_ptr<face> iface = make_shared<intFace>();

      int ic1 = f2c(ff,0);
      // Find local face ID of global face within first element [on left]
      vector<int> cellFaces;
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
        info.isBnd = 0;
        ic1 = eleMap[ic1];
        ic2 = eleMap[ic2];
        if (ic1<0 || ic2<0) {
          FatalError("Unblanking an internal face, but an ele remains blanked!");
        }
        iface->initialize(eles[ic1],eles[ic2],ff,fid1,info,params);
        iface->setupFace();
      }

      // Find the next-lowest index for insertion into vector
      int ind = 0;
      for (int f2=ff-1; f2>=0; f2--) {
        ind = faceMap[f2];
        if ( ind>0 && (faceType[f2] == INTERNAL || faceType[f2] == BOUNDARY) && faces[ind]->ID < ff) {
          ind++;
          break;
        }
      }
      ind = max(ind,0);
      faces.insert(faces.begin()+ind,1,iface);
      faceMap[ff] = ind;
      currFaceType[ff] = INTERNAL;
      nIntFaces++;
      for (int i=ind+1; i<faces.size(); i++)
        faceMap[faces[i]->ID] = i;
    }
    else if (faceType[ff] == BOUNDARY)
    {
      int ind = std::distance(bndFaces.begin(), std::find(bndFaces.begin(),bndFaces.end(),ff));

      if (bcType[ind] == OVERSET) {
        // Boundary face is actually an overset-boundary face
        ubOFaces.insert(ff);
      }
      else {
        // Just a normal boundary face

        shared_ptr<face> bface = make_shared<boundFace>();

        int ic = f2c(ff,0);
        // Find local face ID of global face within element
        vector<int> cellFaces;
        cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
        int fid1 = findFirst(cellFaces,ff);
        if (f2c(ff,1) != -1) {
          FatalError("Boundary face has a right element assigned.");
        }else{
          struct faceInfo info;
          info.bcType = bcType[ind];
          info.isBnd = 1;
          ic = eleMap[ic];
          if (ic<0) FatalError("Unblanking a boundary face, but ele remains blanked!");
          shared_ptr<ele> nullEle;
          bface->initialize(eles[ic],nullEle,ff,fid1,info,params);
          bface->setupFace();
        }

        // Find the next-lowest index for insertion into vector
        int ind = 0;
        for (int f2=ff-1; f2>=0; f2--) {
          ind = faceMap[f2];
          if ( ind>0 && (currFaceType[f2] == INTERNAL || currFaceType[f2] == BOUNDARY) && faces[ind]->ID < ff) {
            ind++;
            break;
          }
        }
        ind = max(ind,0);
        faces.insert(faces.begin()+ind,1,bface);
        faceMap[ff] = ind;
        currFaceType[ff] = BOUNDARY;
        nBndFaces++;
        for (int i=ind+1; i<faces.size(); i++)
          faceMap[faces[i]->ID] = i;
      }
    }
  }

#ifndef _NO_MPI
  for (auto &ff:ubMFaces) {
    int ind = std::distance(mpiFaces.begin(), std::find(mpiFaces.begin(),mpiFaces.end(),ff));

    shared_ptr<mpiFace> mface = make_shared<mpiFace>();

    int ic = f2c(ff,0);
    // Find local face ID of global face within element
    int fid1;
    vector<int> cellFaces;
    if (nDims == 2) {
      cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
      fid1 = findFirst(cellFaces,ff);
    } else {
      fid1 = mpiLocF[ind];
    }

    if (f2c(ff,1) != -1) {
      FatalError("MPI face has a right element assigned.");
    } else {
      int relRot = 0;
      if (nDims == 3) {
        // Find the relative orientation (rotation) between left & right faces
        relRot = compareOrientationMPI(ic,fid1,gIC_R[ind],mpiLocF_R[ind],mpiPeriodic[ind]);
      }
      struct faceInfo info;
      info.IDR = faceID_R[ind];
      info.relRot = relRot;
      info.procL = gridRank;
      info.procR = procR[ind];
      info.isMPI = 1;
      info.gridComm = gridComm;  // Note that this is equivalent to MPI_COMM_WORLD if non-overset (ngrids = 1)
      ic = eleMap[ic];
      if (ic<0) FatalError("Unblanking an MPI face, but ele remains blanked!");
      shared_ptr<ele> nullEle;
      mface->initialize(eles[ic],nullEle,ff,fid1,info,params);
      mface->setupFace();
      ubMFaces.insert(ff);
    }
    // Find the next-lowest index for insertion into vector
    ind = 0;
    for (int i=0; i<mFaces.size(); i++) {
      if (mFaces[i]->ID > ff) {
        ind = i-1;
        break;
      }
    }
    ind = max(ind,0);
    mFaces.insert(mFaces.begin()+ind,1,mface);
    faceMap[ff] = ind;
    currFaceType[ff] = MPI_FACE;
    for (int i=ind+1; i<mFaces.size(); i++) {
      faceMap[mFaces[i]->ID] = i;
    }
  }
#endif

  /* --- Setup Newly-Unblanked Overset Faces --- */

  for (auto &ff:ubOFaces) {
    if (ff<0) continue;

    shared_ptr<overFace> oface = make_shared<overFace>();

    int ic = f2c(ff,0);
    if (ic == -1 || iblankCell[f2c(ff,0)] == HOLE) {
      if (f2c(ff,1) == -1 || iblankCell[f2c(ff,1)] == HOLE) {
        // This happens when a fringe face is ALSO an MPI-boundary face
        // Since the other processor has the non-blanked cell, just ignore the face here
        // But, this should never happen during unblanking... right?
        //continue;
        FatalError("Something went wrong in determining MPI/Overset face unblanking.");
      }
      ic = f2c(ff,1);
    }

    // Find local face ID of global face within first element [on left]
    vector<int> cellFaces;
    cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
    int fid = findFirst(cellFaces,ff);

    struct faceInfo info;
    ic = eleMap[ic];
    if (ic<0) FatalError("Unblanking an Overset face, but ele remains blanked!");
    shared_ptr<ele> nullEle;
    oface->initialize(eles[ic],nullEle,ff,fid,info,params);
    oface->setupFace();

    // Find the next-lowest index for insertion into vector
    int ind = 0;
    while (ind+1<oFaces.size() && oFaces[ind+1]->ID < ff) {
      ind++;
    }
    oFaces.insert(oFaces.begin()+ind,1,oface);
    faceMap[ff] = ind;
    currFaceType[ff] = OVER_FACE;
    for (int i=ind+1; i<oFaces.size(); i++)
      faceMap[oFaces[i]->ID] = i;

    // Add this face to list of overFaces (while keeping list sorted)
    overFaces.insert(ff);
  }
  overFaces.erase(-1);

  // Finish the setup of all unblanked MPI faces (need L/R ranks to be ready)
  for (auto &mface:mFaces)
    if (ubMFaces.count(mface->ID))
      mface->finishRightSetup();
}

void geo::removeFaces(vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces,
                      set<int> &blankIFaces, set<int> &blankMFaces, set<int> &blankOFaces)
{
  // NOTE: faceType refers to the original face type as read from mesh file
  //       (Before overset connectivity processing)

  for (auto &ff:blankIFaces) {
    if (ff<0) continue;

    int ind = faceMap[ff];
    int fType = currFaceType[ff];
    if (ind<0) FatalError("invalid blankIFace!");

    if (fType == INTERNAL || fType == BOUNDARY) {
      faces.erase(faces.begin()+ind,faces.begin()+ind+1);

      faceMap[ff] = -1;
      currFaceType[ff] = HOLE_FACE;
      for (int f2=ff+1; f2<nFaces; f2++) {
        if (currFaceType[f2] == INTERNAL || currFaceType[f2] == BOUNDARY) {
          faceMap[f2]--;
        }
      }

      if (fType == INTERNAL)
        nIntFaces--;
      else
        nBndFaces--;
    }
  }

  for (auto &ff:blankMFaces) {
    if (ff<0) continue;
    int ind = faceMap[ff];
    if (ind<0) FatalError("Invalid blankMFace!");

    mFaces.erase(mFaces.begin()+ind,mFaces.begin()+ind+1);

    faceMap[ff] = -1;
    currFaceType[ff] = HOLE_FACE;
    for (int f2=ff+1; f2<nFaces; f2++) {
      if (currFaceType[f2] == MPI_FACE) {
        faceMap[f2]--;
      }
    }
  }

  for (auto &ff:blankOFaces) {
    if (ff<0) continue;
    int ind = faceMap[ff];
    if (ind<0) continue;

    oFaces.erase(oFaces.begin()+ind,oFaces.begin()+ind+1);

    // Update the map
    faceMap[ff] = -1;
    currFaceType[ff] = HOLE_FACE;
    for (int i=0; i<oFaces.size(); i++) faceMap[oFaces[i]->ID] = i;

    overFaces.erase(ff);
  }
}
