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
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#include "metis.h"
#endif

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
  for (int ib=0; ib<nBounds; ib++) {
    if (bcList[ib] == OVERSET) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iover.push_back(bndPts(ib,iv));
      }
    }
    else if (bcList[ib] == SLIP_WALL || bcList[ib] == ADIABATIC_NOSLIP || bcList[ib] == ISOTHERMAL_NOSLIP) {
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        iwall.push_back(bndPts(ib,iv));
      }
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

  // Create a new communicator to communicate across grids
  MPI_Comm_split(MPI_COMM_WORLD, gridRank, gridID, &interComm);
  MPI_Comm_split(MPI_COMM_WORLD, gridID, params->rank, &gridComm);
#endif
}

void geo::updateOversetConnectivity(void)
{
#ifndef _NO_MPI
  // Pre-process the grids
  tg->profile();

  // Have TIOGA perform the nodal overset connectivity (set nodal iblanks)
  tg->performConnectivity();

  // Only needed for debugging, really
  if (params->writeIBLANK)
    writeOversetConnectivity();

  // Now use new nodal iblanks to set cell and face iblanks
  setCellFaceIblanks();
#endif
}

void geo::setCellFaceIblanks()
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

  // First, blank all cells which contain a hole node

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

  // Only needed for moving grids: List of current hole faces
  holeFaces.clear();
  fringeFaces;
  blankFaces.clear();
  blankOFaces.clear();
  unblankFaces.clear();
  unblankOFaces.clear();
  for (int ff=0; ff<nFaces; ff++) {
    if (iblankFace[ff] == HOLE)
      holeFaces.insert(ff);
    else if (iblankFace[ff] == FRINGE)
      fringeFaces.insert(ff);
  }

  iblankFace.assign(nFaces,NORMAL);

  // Next, get the new overset faces & set all hole faces

  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic] == HOLE) {
      for (int j=0; j<c2nf[ic]; j++) {
        int ff = c2f(ic,j);
        if (iblankFace[ff] == NORMAL) {
          // If not set yet, assume fringe (overset), but set to hole if any hole nodes
          iblankFace[ff] = FRINGE;
          for (int k=0; k<f2nv[ff]; k++) {
            int iv = f2v(ff,k);
            if (iblank[iv] == HOLE) {
              iblankFace[ff] = HOLE;
              break;
            }
          }
        }
      }
    }
  }

  // Figure out which faces need to be removed from where
  // ('Normal' face going to either 'fringe' or 'hole', or 'fringe' faces going to 'hole')
  for (int ff=0; ff<nFaces; ff++) {
    if (!holeFaces.count(ff)) {
      // Not a hole face, so it exists somewhere
      if (!fringeFaces.count(ff)) {
        // Not a fringe face either; 'normal' face (int, bound, or mpi)
        if (iblankFace[ff] != NORMAL)
          blankFaces.insert(ff);
      }
      else {
        // Fringe face; if no longer fringe, must be removed so mark as such
        if (iblankFace[ff] != FRINGE)
          blankOFaces.insert(ff);
      }
    }
  }

  // Only needed for moving grids: Get faces which  must be 'un-blanked'
  for (auto &ff:holeFaces) {
    if (iblankFace[ff] == NORMAL) {
      // Face will be created as an int, bound, or mpi face
      unblankFaces.insert(ff);
    }
    else if (iblankFace[ff] == FRINGE) {
      // Face will be created as an overset face
      unblankOFaces.insert(ff);
    }
  }
}

void geo::writeOversetConnectivity(void)
{
#ifndef _NO_MPI
  // Write out only the mesh with IBLANK info (no solution data)
  tg->writeData(0,NULL,0);
#endif
}

void geo::matchOversetDonors(vector<ele> &eles, vector<superMesh> &donors)
{
#ifndef _NO_MPI

#endif
}

void geo::setupUnblankElesFaces(vector<ele> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces)
{
  /* --- Remove Newly-Blanked Elements --- */

  for (auto &ic:blankCells) {
    int ind = eleMap[ic];
    eles.erase(eles.begin()+ind,eles.begin()+ind+1);
    eleMap[ic] = -1;

    // Update the map
    for (int k=ic+1; k<nEles; k++)
      eleMap[k]--;
  }

  /* --- Setup & Insert Unblanked Elements --- */

  for (auto &ic:unblankCells) {
    // Find the next-lowest index
    int ind = eleMap[ic];
    int j = 0;
    while (ind < 0) {
      ind = eleMap[ic-j];
      j++;
    }
    ind++;

    ele e;
    e.ID = ic;
    if (nprocPerGrid>1)
      e.IDg = ic2icg[ic];
    else
      e.IDg = ic;
    e.eType = ctype[ic];
    e.nNodes = c2nv[ic];
    if (nDims == 2)
      e.nMpts = 4;
    else
      e.nMpts = 8;

    // Shape [mesh] nodes
    e.nodeID.resize(c2nv[ic]);
    e.nodes.resize(c2nv[ic]);
    for (int iv=0; iv<c2nv[ic]; iv++) {
      e.nodeID[iv] = c2v(ic,iv);
      e.nodes[iv] = point(xv[c2v(ic,iv)]);
    }

    // Global face IDs for internal & boundary faces
    e.faceID.resize(c2nf[ic]);
    e.bndFace.resize(c2nf[ic]);
    for (int k=0; k<c2nf[ic]; k++) {
      e.bndFace[k] = c2b(ic,k);
      e.faceID[k] = c2f(ic,k);
    }

    eles.insert(eles.begin()+ind,1,e);

    // Update the map
    eleMap[ic] = ind;
    for (int k=ic+1; k<nEles; k++)
      eleMap[k]++;
  }

  /* --- Remove Newly-Blanked Faces --- */

  for (auto &ff:blankFaces) {
    int ind = faceMap[ff];
    int fType = faceType[ff];

    if (fType == INTERNAL || fType == BOUNDARY) {
      faces.erase(faces.begin()+ind,faces.begin()+ind+1);

      faceMap[ff] = -1;
      for (int f2=ff+1; f2<nFaces; f2++) {
        if (!fringeFaces.count(f2) && (faceType[f2] == INTERNAL || faceType[f2] == BOUNDARY)) {
          faceMap[f2]--;
        }
      }

      if (fType == INTERNAL)
        nIntFaces--;
      else
        nBndFaces--;
    }
    else if (fType == MPI_FACE) {
      mFaces.erase(mFaces.begin()+ind,mFaces.begin()+ind+1);

      faceMap[ff] = -1;
      for (int f2=ff+1; f2<nFaces; f2++) {
        if (!fringeFaces.count(f2) && faceType[f2] == fType) {
          faceMap[f2]--;
        }
      }
    }
    else {
      FatalError("Face does not have a proper type assigned!");
    }
  }

  for (auto &ff:blankOFaces) {
    int ind = faceMap[ff];
    oFaces.erase(oFaces.begin()+ind,oFaces.begin()+ind+1);

    // Update the map
    faceMap[ff] = -1;
    for (auto &f2:fringeFaces) {
      if (f2 > ff)
        faceMap[f2]--;
    }
  }

  /* --- Create Newly-Unblanked Faces --- */

  for (auto &ff:unblankFaces) {
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
        ic1 = eleMap[ic1];
        ic2 = eleMap[ic2];
        iface->initialize(&eles[ic1],&eles[ic2],ff,fid1,info,params);
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
      faces.insert(faces.begin()+ind,1,iface);
    }
    else if (faceType[ff] == BOUNDARY)
    {
      int ind = std::distance(bndFaces.begin(), std::find(bndFaces.begin(),bndFaces.end(),ff));

      if (bcType[ind] == OVERSET) {
        // Boundary face is actually an overset-boundary face
        unblankOFaces.insert(ff);
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
          bface->initialize(&eles[ic],NULL,ff,fid1,info,params);
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
        faces.insert(faces.begin()+ind,1,bface);
      }
    }
    else if (faceType[ff] == MPI_FACE)
    {
      int ind = std::distance(mpiFaces.begin(), std::find(mpiFaces.begin(),mpiFaces.end(),ff));

      shared_ptr<mpiFace> mface = make_shared<mpiFace>();

      int ic = f2c(ff,0);
      // Find local face ID of global face within element
      int fid1;
      vector<int> cellFaces;
      if (nDims == 2) {
        cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
        fid1 = findFirst(cellFaces,ff);
      }
      else {
        fid1 = mpiLocF[ind];
      }
      if (f2c(ff,1) != -1) {
        FatalError("MPI face has a right element assigned.");
      }else{
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
        mface->initialize(&eles[ic],NULL,ff,fid1,info,params);
      }
    }
  }

  /* --- Setup Newly-Unblanked Overset Faces --- */

  for (auto &ff:unblankOFaces) {
    shared_ptr<overFace> oface = make_shared<overFace>();

    int ic = f2c(ff,0);
    if (ic == -1 || iblankCell[f2c(ff,0)] == HOLE) {
      if (f2c(ff,1) == -1 || iblankCell[f2c(ff,1)] == HOLE) {
        // This happens when a fringe face is ALSO an MPI-boundary face
        // Since the other processor has the non-blanked cell, just ignore the face here
        //ff = -1; // to remove from vector later
        continue;
      }
      ic = f2c(ff,1);
    }

    // Find local face ID of global face within first element [on left]
    vector<int> cellFaces;
    cellFaces.assign(c2f[ic],c2f[ic]+c2nf[ic]);
    int fid = findFirst(cellFaces,ff);

    struct faceInfo info;
    ic = eleMap[ic];
    oface->initialize(&eles[ic],NULL,ff,fid,info,params);
  }
}
