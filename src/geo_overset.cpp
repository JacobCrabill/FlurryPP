/*!
 * \file geo_overset.cpp
 * \brief Overset-related methods for the geo class
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
#include <unordered_set>
#include <sstream>

#include "ele.hpp"
#include "face.hpp"
#include "overFace.hpp"
#include "overComm.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "ADT.h"
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

  vector<int> nProcsGrid(nGrids);
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

  gridIdList.resize(nproc);
  MPI_Allgather(&gridID,1,MPI_INT,gridIdList.data(),1,MPI_INT,MPI_COMM_WORLD);
#endif
}

void geo::setupOverset3D(void)
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
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        nodeType[bndPts(ib,iv)] = BOUNDARY_NODE;
      }
    }
  }

  int nwall = iwall.size();
  int nover = iover.size();
  int ntypes = 1;                  //! Number of element types in grid block
  nodesPerCell = new int[ntypes];
  nodesPerCell[0] = nNodesPerCell; //! Number of nodes per element for each element type (but only one type so far)
  iblank.resize(nVerts);
  iblankCell.resize(nEles);
  iblankFace.resize(nFaces);

  // Need an int**, even if only have one element type
  conn[0] = c2v.getData();

  tg->registerGridData(gridID,nVerts,xv.getData(),iblank.data(),nwall,nover,iwall.data(),
                       iover.data(),ntypes,nodesPerCell,&nEles,&conn[0]);

  tg->set_cell_iblank(iblankCell.data());

  // Get a list of all cells which have an overset-boundary face (for later use with blanking)
  for (int i=0; i<nBndFaces; i++) {
    if (bcType[i]==OVERSET) {
      overCells.insert(f2c(bndFaces[i],0));
    }
  }

  // Pre-process the grids
  tg->profile();

  // Have TIOGA perform the nodal overset connectivity (set nodal iblanks)
  tg->performConnectivity();

  // Now use new nodal iblanks to set cell and face iblanks
  setIblankEles(iblank,iblankCell);
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

  eleBBox.setup(nEles,nDims*2);
  for (int i=0; i<nEles; i++) {
    matrix<double> epts;
    for (int j=0; j<c2nv[i]; j++) {
      epts.insertRow(xv[c2v(i,j)],INSERT_AT_END,nDims);
    }
    getBoundingBox(epts.getData(),c2nv[i],nDims,eleBBox[i]);
  }

  adt = make_shared<ADT>();
  OComm->adt = adt;

  adt->buildADT(nDims*2,nEles,eleBBox.getData());

  // Now use new nodal iblanks to set cell and face iblanks
  setIblankEles(iblank,iblankCell);

  fringeCells.clear();
  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic] == FRINGE)
      fringeCells.insert(ic);
  }
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
      for (int iv=0; iv<nBndPts[ib]; iv++) {
        nodeType[bndPts(ib,iv)] = BOUNDARY_NODE;
      }
    }
  }

  unordered_set<int> overNodes;
  overFaceNodes.setup(0,0);
  wallFaceNodes.setup(0,0);
  for (int bf=0; bf<nBndFaces; bf++) {
    if (bcType[bf] == SLIP_WALL || bcType[bf] == ADIABATIC_NOSLIP || bcType[bf] == ISOTHERMAL_NOSLIP) {
      int ff = bndFaces[bf];
      wallFaceNodes.insertRow({f2v(ff,0),f2v(ff,1)});
    }
    else if (bcType[bf] == OVERSET) {
      int ff = bndFaces[bf];
      overFaces.insert(ff);
      if (params->oversetMethod != 0) {
        overNodes.insert(f2v(ff,0));
        overNodes.insert(f2v(ff,1));
      } else {
        overFaceNodes.insertRow({f2v(ff,0),f2v(ff,1)});
      }
    }
  }

  if (params->oversetMethod != 0) {
    unordered_set<int> newOverFaces;
    // Need more overlap for field interp method; move back boundary by one layer
    for (auto &ff:overFaces) {
      int ic = f2c(ff,0);
      if (params->oversetMethod == 2)
        fringeCells.insert(ic);
      for (int j=0; j<c2nf[ic]; j++) {
        int ff2 = c2f(ic,j);
        int iv1 = f2v(ff2,0);
        int iv2 = f2v(ff2,1);
        if (!overNodes.count(iv1) && !overNodes.count(iv2)) {
          if (params->oversetMethod == 2)
            newOverFaces.insert(ff2);
          overFaceNodes.insertRow({iv1,iv2});
        }
      }
    }
    if (params->oversetMethod == 2) {
      overFaces.clear();
      overFaces = newOverFaces;
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
  } else {
    for (int i=0; i<nEles; i++) {
      matrix<double> epts;
      for (int j=0; j<c2nv[i]; j++) {
        epts.insertRow(xv[c2v(i,j)],INSERT_AT_END,nDims);
      }
      getBoundingBox(epts.getData(),c2nv[i],nDims,eleBBox[i]);
    }
    adt->buildADT(nDims*2,nEles,eleBBox.getData());
  }
#endif
}

void geo::setIterIblanks(void)
{
#ifndef _NO_MPI
  if (params->iter==params->initIter) {
    holeCells.clear();
    for (int ic=0; ic<nEles; ic++)
      if (iblankCell[ic] == HOLE)
        holeCells.insert(ic);
  }

  /* ---- Get iblank data for end of iteration ---- */

  moveMesh(1.);

  vector<int> iblankVert1;
  if (nDims == 2) {
    OComm->setIblanks2D(xv,overFaceNodes,wallFaceNodes,iblankVert1);
  } else {
    // Pre-process the grids, then have TIOGA set nodal iblank values
    tg->profile();
    tg->performConnectivity();
    iblankVert1 = iblank;
  }

  //! TODO: in 3D - need better hole blanking...
  if (nDims == 3) {
    for (int iv = 0; iv < nVerts; iv++)
      if (nodeType[iv] == OVERSET_NODE && iblankVert1[iv] != HOLE)
        iblankVert1[iv] = NORMAL;
  }

  /* ---- Get iblank data for beginning of iteration ---- */

  moveMesh(0.);

  if (nDims == 2) {
    OComm->setIblanks2D(xv,overFaceNodes,wallFaceNodes,iblank);
  } else {
    // Pre-process the grids, then have TIOGA set nodal iblank values
    tg->profile();
    tg->performConnectivity();
  }

  //! TODO: in 3D - need better hole blanking...
  if (nDims == 3) {
    for (int iv = 0; iv < nVerts; iv++) {
      if (nodeType[iv] == OVERSET_NODE && iblank[iv] != HOLE)
        iblank[iv] = NORMAL;
    }
  }

  // Take the union of the hole and normal regions, leaving the intersection of
  // the fringe regions
  for (int iv=0; iv<nVerts; iv++) {
    if (iblankVert1[iv] == HOLE || iblank[iv] == HOLE)
      iblank[iv] = HOLE;

    if (iblankVert1[iv] == NORMAL || iblank[iv] == NORMAL)
      iblank[iv] = NORMAL;
  }

  setIblankEles(iblank,iblankCell);

  unordered_set<int> holeCells_tmp;
  fringeCells.clear();
  blankCells.clear();
  unblankCells.clear();
  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic] == HOLE) {
      holeCells_tmp.insert(ic);
      if (!holeCells.count(ic))
        blankCells.insert(ic);
    } else {
      if (iblankCell[ic] == FRINGE)
        fringeCells.insert(ic);
      if (holeCells.count(ic))
        unblankCells.insert(ic);
    }
  }
  holeCells = holeCells_tmp;
#endif
}

void geo::setIblankEles(vector<int> &iblankVert, vector<int> &iblankEle)
{
#ifndef _NO_MPI
  // Use the TIOGA-supplied nodal iblank values, set iblank values for all cells and faces

  // Only needed for moving grids: List of current hole cells
  iblankEle.assign(nEles,NORMAL);

  /* --- Move on to cell iblank values --- */

  // Simply blank all cells which contain a hole node
  for (int ic=0; ic<nEles; ic++) {
    for (int j=0; j<c2nv[ic]; j++) {
      int iv = c2v(ic,j);
      if (iblankVert[iv] == HOLE) {
        iblankEle[ic] = HOLE;
        break;
      }
    }
  }

  vector<int> mpiFringeFaces;
  vector<int> iblankEle1;
  if (params->oversetMethod != 2)
    iblankEle1 = iblankEle;
  else
    iblankEle1.assign(nEles,NORMAL);

  if (params->oversetMethod != 2) {
    // For either of the Artificial Boundary methods
    for (int ic=0; ic<nEles; ic++) {
      if (iblankEle[ic] == NORMAL) {
        int nfringe = 0;
        for (int j=0; j<c2nv[ic]; j++) {
          if (iblankVert[c2v(ic,j)] == FRINGE) {
            nfringe++;
          }
        }
        if (nfringe == c2nv[ic]) {
          iblankEle1[ic] = FRINGE;
          iblankEle[ic]  = FRINGE;
        }
      }
    }

    if (params->oversetMethod == 1) {
      // Extra overlap needed - Bring back outermost 'fringe' hole-cut layer as 'normal' cells
      for (int ic=0; ic<nEles; ic++) {
        if (iblankEle[ic] == NORMAL) {
          for (int j=0; j<c2nf[ic]; j++) {
            int ic2 = c2c(ic,j);
            if (ic2 > -1) {
              if (iblankEle1[ic2] == FRINGE) {
                iblankEle1[ic2] = NORMAL;
              }
            } else {
              // MPI Boundary
              int F = findFirst(mpiFaces,c2f(ic,j));
              if (F > -1) {
                mpiFringeFaces.push_back(mpiFaces[F]);
              }
            }
          }
        }
      }
    }
    else if (params->oversetMethod == 0 && nDims == 3) {
      for (int ic=0; ic<nEles; ic++) {
        if (iblankEle[ic] != NORMAL) {
          for (int j=0; j<c2nf[ic]; j++) {
            int ic2 = c2c(ic,j);
            if (ic2 > -1) {
              if (iblankEle1[ic2] == NORMAL) {
                iblankEle1[ic2] = FRINGE;
              }
            } else {
              // MPI Boundary
              int F = findFirst(mpiFaces,c2f(ic,j));
              if (F > -1) {
                mpiFringeFaces.push_back(mpiFaces[F]);
              }
            }
          }
        }
      }
    }
  }

  else if (params->oversetMethod == 2) {
    // Tag the innermost layer of 'normal' cells as 'fringe' cells for unblank method
    for (int ic=0; ic<nEles; ic++) {
      if (iblankEle[ic] == NORMAL) {
        int nfringe = 0;
        for (int j=0; j<c2nv[ic]; j++) {
          if (iblankVert[c2v(ic,j)] == FRINGE) {
            nfringe++;
          }
        }
        if (nfringe == c2nv[ic]) {
          iblankEle[ic] = FRINGE;
        }
      }
    }

    // Step back by one layer of cells

    for (int ic=0; ic<nEles; ic++) {
      if (iblankEle[ic] == NORMAL) {
        for (int j=0; j<c2nf[ic]; j++) {
          int ic2 = c2c(ic,j);
          if (ic2 > -1) {
            if (iblankEle[ic2] == FRINGE) {
              iblankEle1[ic2] = FRINGE;
            }
          } else {
            // MPI Boundary
            int F = findFirst(mpiFaces,c2f(ic,j));
            if (F > -1) {
              mpiFringeFaces.push_back(mpiFaces[F]);
            }
          }
        }
      }
    }
  }

  if (params->oversetMethod > 0 || nDims == 3) {
    // Enforce consistency across MPI boundaries for fringe->hole conversion
    vector<int> nFringe_proc(nProcGrid);
    int nFringe = mpiFringeFaces.size();

    MPI_Allgather(&nFringe,1,MPI_INT,nFringe_proc.data(),1,MPI_INT,gridComm);

    int maxNFringe = getMax(nFringe_proc);
    matrix<int> mpiFringeFaces_proc(nProcGrid,maxNFringe);

    vector<int> recvCnts(nProcGrid);
    vector<int> recvDisp(nProcGrid);
    for (int i=0; i<nProcGrid; i++) {
      recvCnts[i] = nFringe_proc[i];
      recvDisp[i] = i*maxNFringe;
    }
    MPI_Allgatherv(mpiFringeFaces.data(),mpiFringeFaces.size(),MPI_INT,mpiFringeFaces_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

    if (params->oversetMethod == 1) {
      for (int F=0; F<nMpiFaces; F++) {
        int ff = mpiFaces[F];
        int fr = faceID_R[F];
        if (iblankEle1[f2c(ff,0)] == FRINGE) {
          int p = procR[F];
          for (int i=0; i<recvCnts[p]; i++) {
            int fr2 = mpiFringeFaces_proc(p,i);
            if (fr == fr2) {
              iblankEle1[f2c(ff,0)] = NORMAL;
              break;
            }
          }
        }
      }
    }
    else if (params->oversetMethod == 0 && nDims == 3) {
            for (int F=0; F<nMpiFaces; F++) {
        int ff = mpiFaces[F];
        int fr = faceID_R[F];
        if (iblankEle1[f2c(ff,0)] == NORMAL) {
          int p = procR[F];
          for (int i=0; i<recvCnts[p]; i++) {
            int fr2 = mpiFringeFaces_proc(p,i);
            if (fr == fr2) {
              iblankEle1[f2c(ff,0)] = FRINGE;
              break;
            }
          }
        }
      }
    }

    else if (params->oversetMethod == 2) {
      for (int F=0; F<nMpiFaces; F++) {
        int ff = mpiFaces[F];
        int fr = faceID_R[F];
        if (iblankEle[f2c(ff,0)] == FRINGE) {
          int p = procR[F];
          for (int i=0; i<recvCnts[p]; i++) {
            int fr2 = mpiFringeFaces_proc(p,i);
            if (fr == fr2) {
              iblankEle1[f2c(ff,0)] = FRINGE;
              break;
            }
          }
        }
      }

      for (int ic=0; ic<nEles; ic++) {
        if (iblankEle1[ic] == NORMAL && iblankEle[ic] == FRINGE) {
          iblankEle[ic] = HOLE;
        }
      }

      for (int i=0; i<bndFaces.size(); i++) {
        if (bcType[i] == OVERSET) {
          int ic = f2c(bndFaces[i],0);
          iblankEle[ic] = FRINGE;
        }
      }
    }
  }

  //! ---- QUICK HACK FOR 3D SPHERE CASE ----
  if (params->oversetMethod == 0 && params->motion==5 && nDims == 3) {
    for (int ic = 0; ic < nEles; ic++) {
      if (iblankEle[ic] == NORMAL && iblankEle1[ic] != NORMAL)
        iblankEle[ic] = iblankEle1[ic];
      else if (iblankEle[ic] == FRINGE && iblankEle1[ic] == HOLE)
        iblankEle[ic] = HOLE;
    }
    iblankEle1 = iblankEle;

    mpiFringeFaces.resize(0);
    for (int ic=0; ic<nEles; ic++) {
      if (iblankEle[ic] != NORMAL) {
        for (int j=0; j<c2nf[ic]; j++) {
          int ic2 = c2c(ic,j);
          if (ic2 > -1) {
            if (iblankEle1[ic2] == NORMAL) {
              iblankEle1[ic2] = FRINGE;
            }
          } else {
            // MPI Boundary
            int F = findFirst(mpiFaces,c2f(ic,j));
            if (F > -1) {
              mpiFringeFaces.push_back(mpiFaces[F]);
            }
          }
        }
      }
    }

    // Enforce consistency across MPI boundaries for fringe->hole conversion
    vector<int> nFringe_proc(nProcGrid);
    int nFringe = mpiFringeFaces.size();
    MPI_Allgather(&nFringe,1,MPI_INT,nFringe_proc.data(),1,MPI_INT,gridComm);

    int maxNFringe = getMax(nFringe_proc);
    matrix<int> mpiFringeFaces_proc(nProcGrid,maxNFringe);

    vector<int> recvCnts(nProcGrid);
    vector<int> recvDisp(nProcGrid);
    for (int i=0; i<nProcGrid; i++) {
      recvCnts[i] = nFringe_proc[i];
      recvDisp[i] = i*maxNFringe;
    }
    MPI_Allgatherv(mpiFringeFaces.data(),mpiFringeFaces.size(),MPI_INT,mpiFringeFaces_proc.getData(),recvCnts.data(),recvDisp.data(),MPI_INT,gridComm);

    for (int F=0; F<nMpiFaces; F++) {
      int ff = mpiFaces[F];
      int fr = faceID_R[F];
      if (iblankEle1[f2c(ff,0)] == NORMAL) {
        int p = procR[F];
        for (int i=0; i<recvCnts[p]; i++) {
          int fr2 = mpiFringeFaces_proc(p,i);
          if (fr == fr2) {
            iblankEle1[f2c(ff,0)] = FRINGE;
            break;
          }
        }
      }
    }
  }
  //! ---- END HACK ----

  // Final blanking update for both AB methods
  if (params->oversetMethod != 2) {
    for (int ic=0; ic<nEles; ic++) {
      if (iblankEle1[ic] == NORMAL)
        iblankEle[ic] = NORMAL;
      else
        iblankEle[ic] = HOLE;
    }
  }

  ///! EVEN BIGGER HACK
  if (nDims == 3 && params->motion==5 && gridID == 0) iblankEle.assign(nEles,NORMAL);
#endif
}

void geo::updateBlankingTioga(void)
{
#ifndef _NO_MPI
  if (nDims == 3) {
    // Pre-process the grids
    tg->profile();

    // Have TIOGA perform the nodal overset connectivity (set nodal iblanks)
    tg->performConnectivity();

    // Now use new nodal iblanks to set cell and face iblanks
    setIblankEles(iblank,iblankCell);
  }
#endif
}

void geo::setFaceIblanks(void)
{
  overFaces.clear();
  iblankFace.assign(nFaces,NORMAL);

  for (int ic=0; ic<nEles; ic++) {
    if (iblankCell[ic]!=HOLE) continue;

    for (int j=0; j<c2nf[ic]; j++) {
      int ff = c2f(ic,j);
      if (c2c(ic,j)>=0) {
        // Internal face
        if (iblankCell[c2c(ic,j)] == HOLE) {
          iblankFace[ff] = HOLE;
        } else {
          iblankFace[ff] = FRINGE;
          overFaces.insert(ff);
        }
      } else {
        // Boundary or MPI face
        iblankFace[ff] = HOLE;
      }
    }
  }

#ifndef _NO_MPI
  // Get the number of mpiFaces on each processor (for later communication)
  vector<int> nMpiFaces_proc(nProcGrid);
  MPI_Allgather(&nMpiFaces,1,MPI_INT,nMpiFaces_proc.data(),1,MPI_INT,gridComm);
  int maxNMpiFaces = getMax(nMpiFaces_proc);

  vector<int> mpiIblank(nMpiFaces);
  matrix<int> mpiIblank_proc(nProcGrid,maxNMpiFaces);
  matrix<int> mpiFid_proc(nProcGrid,maxNMpiFaces);

  vector<int> recvCnts(nProcGrid);
  vector<int> recvDisp(nProcGrid);
  for (int i=0; i<nProcGrid; i++) {
    recvCnts[i] = nMpiFaces_proc[i];
    recvDisp[i] = i*maxNMpiFaces;
  }

  for (int i = 0; i < nMpiFaces; i++)
    mpiIblank[i] = iblankCell[f2c(mpiFaces[i],0)];

  // Get iblank data for all mpi faces
  MPI_Allgatherv(mpiFaces.data(), nMpiFaces, MPI_INT, mpiFid_proc.getData(), recvCnts.data(), recvDisp.data(), MPI_INT, gridComm);
  MPI_Allgatherv(mpiIblank.data(), nMpiFaces, MPI_INT, mpiIblank_proc.getData(), recvCnts.data(), recvDisp.data(), MPI_INT, gridComm);

  for (int F = 0; F < nMpiFaces; F++) {
    int ff = mpiFaces[F];
    int p = procR[F];
    for (int i = 0; i < nMpiFaces_proc[p]; i++) {
      if (mpiFid_proc(p,i) == faceID_R[F]) {
        if (mpiIblank[F] != NORMAL || mpiIblank_proc(p,i) != NORMAL) {
          // Not a normal face; figure out if hole or fringe
          if (mpiIblank[F] == HOLE && mpiIblank_proc(p,i) == HOLE)
            iblankFace[ff] = HOLE;
          else {
            iblankFace[ff] = FRINGE;
            overFaces.insert(ff);
          }
        }
      }
    }
  }

#endif
}

void geo::matchOversetDonors(vector<shared_ptr<ele>> &eles, vector<superMesh> &donors)
{
#ifndef _NO_MPI

#endif
}

void geo::processBlanks(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces, solver *Solver)
{
#ifndef _NO_MPI
  /* --- Check whether anything needs to be done --- */

  int NB = blankCells.size();
  int nBlanks;
  MPI_Allreduce(&NB,&nBlanks,1,MPI_INT,MPI_SUM,gridComm);

  if (nBlanks == 0) return;

  /* --- Set blank/unblank faces for all elements to be blanked --- */

  unordered_set<int> blankIFaces, blankMFaces, blankOFaces, ubOFaces;
  for (auto &ic:blankCells) {
    if (eleMap[ic] < 0) continue; // Ignore already-blanked cells

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
    if ( (iblankCell[ic1]!=HOLE && iblankCell[ic2]==HOLE) ||
         (iblankCell[ic1]==HOLE && iblankCell[ic2]!=HOLE) )
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

  unordered_set<int> ubIFaces, ubMFaces;

  removeEles(eles,blankCells,Solver);
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

void geo::processUnblanks(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces, solver *Solver)
{
#ifndef _NO_MPI
  /* --- Check whether anything needs to be done --- */

  int nUB = unblankCells.size();
  int nUnblanks;
  MPI_Allreduce(&nUB,&nUnblanks,1,MPI_INT,MPI_SUM,gridComm);

  if (nUnblanks == 0) return;

  /* --- Set Unblank/Blank Faces for All Unblank Elements --- */

  unordered_set<int> ubIntFaces, ubMpiFaces, ubOFaces;
  unordered_set<int> blankIFaces, blankMFaces, blankOFaces;
  for (auto &ic:unblankCells) {
    for (int j=0; j<c2nf[ic]; j++) {
      int ic2 = c2c(ic,j);
      int ff2 = c2f(ic,j);
      if (ic2>=0) {
        if (iblankCell[ic2]!=HOLE || unblankCells.count(ic2))
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

  // Figure out whether MPI faces are to be unblanked as MPI or as overset
  // Need list of all possibly-unblanked MPI faces on this rank
  unordered_set<int> allMpiFaces = ubMpiFaces;
  for (auto &ff:overFaces) if (faceType[ff]==MPI_FACE) allMpiFaces.insert(ff);

  vector<int> ubMpi;
  for (auto &ff:allMpiFaces) ubMpi.push_back(ff);

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

  // Now, figure out what faces must be removed due to being replaced by other type
  // For cell unblanking, the only possibility for face blanking is overset faces

  for (auto &ff:ubIntFaces)
    if (overFaces.count(ff))
      blankOFaces.insert(ff);

  for (auto &ff:ubMpiFaces)
    if (overFaces.count(ff))
      blankOFaces.insert(ff);

  insertEles(eles,unblankCells,Solver);
  removeFaces(faces,mFaces,oFaces,blankIFaces,blankMFaces,blankOFaces);
  insertFaces(eles,faces,mFaces,oFaces,ubIntFaces,ubMpiFaces,ubOFaces);

  for (auto &iface:faces) {
    iface->getPointers();
    iface->getPointersRight();
  }
  for (auto &mface:mFaces) mface->getPointers();
  for (auto &oface:oFaces) oface->getPointers();
#endif
}

void geo::removeEles(vector<shared_ptr<ele>> &eles, unordered_set<int> &blankEles, solver *Solver)
{
  /* --- Remove Newly-Blanked Elements --- */

  for (auto &ic:blankEles) {
    if (ic<0) continue;

    int ind = eleMap[ic];
    if (ind<0) continue; //FatalError("Should not have marked a hole cell for blanking!");
    eles.erase(eles.begin()+ind,eles.begin()+ind+1);
    eleMap[ic] = -1;

    Solver->removeElement(ind);

    // Update the map
    for (int k=ic+1; k<nEles; k++)
      if (eleMap[k]>=0)
        eleMap[k]--;

    for (int k = ind; k < eles.size(); k++) {
      eles[k]->sID = k;
    }
  }
}

void geo::insertEles(vector<shared_ptr<ele>> &eles, unordered_set<int> &ubEles, solver *Solver)
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
    e->sID = ind;
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

    Solver->insertElement(ind);

    e->setup(params,Solver,this);

    eles.insert(eles.begin()+ind,1,e);

    // Update the map
    eleMap[ic] = ind;
    for (int k = ic+1; k < nEles; k++)
      if (eleMap[k]>=0)
        eleMap[k]++;

    for (int k = ind+1; k < eles.size(); k++)
      eles[k]->sID = k;
  }

  if (ubEles.size() > 0)
    Solver->updatePosSptsFpts();
}

void geo::insertFaces(vector<shared_ptr<ele>> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mFaces, vector<shared_ptr<overFace>> &oFaces,
                      unordered_set<int> &ubIFaces, unordered_set<int> &ubMFaces, unordered_set<int> &ubOFaces)
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
      for (ind=0; ind<faces.size(); ind++)
        if (faces[ind]->ID > ff)
          break;

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
          if ( ind>0 && (currFaceType[f2] == INTERNAL || currFaceType[f2] == BOUNDARY) && faces[ind]->ID<ff ) {
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
        ind = i;
        break;
      }
    }
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
    for (int i=0; i<oFaces.size(); i++) {
      if (oFaces[i]->ID > ff) {
        ind = i;
        break;
      }
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
                      unordered_set<int>& blankIFaces, unordered_set<int>& blankMFaces, unordered_set<int>& blankOFaces)
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
      for (int i=0; i<faces.size(); i++) faceMap[faces[i]->ID] = i;

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
    for (int i=0; i<mFaces.size(); i++) faceMap[mFaces[i]->ID] = i;
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
