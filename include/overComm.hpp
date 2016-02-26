/*!
 * \file overComm.hpp
 * \brief Header file for overComm class
 *
 * Handles communication of data across multiple MPI-partitioned overset grids
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
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA..
 *
 */
#pragma once

#include "global.hpp"

#include <map>
#include <set>
#include <unordered_set>
#include <vector>

class oper;
class overFace;

#include "ele.hpp"
#include "input.hpp"
#include "operators.hpp"
#include "overFace.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "tioga.h"
#include "ADT.h"
#include "mpi.h"
#endif

#ifndef _NO_MPI
template<typename T>
MPI_Datatype getMpiDatatype(void);
#endif

class overComm
{
public:
  overComm();

  int nGrids;
  int nprocPerGrid;
  int gridID;
  int gridRank;
  int rank;
  int nproc;

  vector<int> gridIdList;  //! Grid ID for each MPI rank

  int nFields;

  input *params;  //! Simulation input parameters

#ifndef _NO_MPI
  MPI_Comm gridComm;

  shared_ptr<tioga> tg;  //! TIOGA object in use for simulation
  shared_ptr<ADT> adt;   //! Alternating Digital Tree for searching
#endif

  /* --- Variables for Exchanging Data at Overset Faces --- */

  vector<int> nPts_rank;           //! Number of fringe points for each rank of current grid
  vector<vector<int>> foundPts;    //! IDs of receptor points from each grid which were found to lie within current grid
  vector<vector<int>> foundRank;   //! gridRank of this process for each found point (for benefit of other processes; probably not needed)
  vector<vector<int>> foundEles;   //! Ele ID which each matched point was found to lie within
  vector<vector<point>> foundLocs; //! Reference location within donor ele of each matched receptor point
  vector<vector<point>> foundNorm; //! For corrected-flux method: Outward unit normal at fringe boundary point

  int nOverPts;                 //! Number of overset (receptor) points on this grid
  matrix<double> overPts;       //! Physical locations of the receptor points on this grid
  matrix<double> overNorm;      //! Outward unit normals at fringe-boundary points on this grid [corrected-flux interp method]
  vector<int> nPtsRecv;         //! Number of points incoming from each grid (across interComm)
  vector<int> nPtsSend;         //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvPts;  //! Point IDs which will be received from each grid (across interComm) (counter to foundPts)

  matrix<double> U_in;          //! Data received from other grid(s)
  vector<matrix<double>> U_out; //! Interpolated data being sent to other grid(s)

  matrix<double> gradU_in;          //! Gradient data received from other grid(s)
  vector<matrix<double>> gradU_out; //! Interpolated gradient data being sent to other grid(s)

  /* --- Variables for Exchanging Data on Unblanked Cells --- */

  vector<int> nCells_rank;              //! Number of unblanked cells for each rank of current grid
  vector<vector<int>> foundCells;       //! IDs of unblanked cells from each grid which were found to overlap this grid
  vector<matrix<int>> foundCellDonors;  //! Donor cells on this grid for each found unblanked cell
  vector<vector<int>> foundCellNDonors; //! Number of donor cells on this grid for each found unblanked cell

  int nUnblanks;          //! Number of unblank cells on this grid
  int nUnblanksTotal;     //! Total number of cells to unblank across domain
  vector<int> ubCells;    //! Cells from this grid which need to be unblanked
  vector<int> nCellsRecv;         //! Number of points incoming from each grid (across interComm)
  vector<int> nCellsSend;         //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvCells;  //! Cell IDs which will be received from each grid (across interComm) (counter to foundPts)

  vector<int> nQptsSend;
  vector<int> nQptsRecv;

  //! Local supermesh of donor elements for each cell needing to be unblanked
  vector<superMesh> donors;

  /* --- Member Functions --- */

  void setup(input *_params, int _nGrids, int _gridID, int _gridRank, int _nprocPerGrid, vector<int> &_gridIdList);

  void setIblanks2D(matrix<double> &xv, matrix<int>& overFaces, matrix<int> &wallFaces, vector<int> &iblank);

  //! Find all elements from eles which overlap with target bounding-box
  unordered_set<int> findCellDonors2D(vector<shared_ptr<ele>> &eles, const vector<double> &targetBox);

  /*!
   * \brief Match up each overset-face flux point to its donor grid and element
   *
   * @param[in] c2ac    : List (for each cell in grid partition) of all cells which share at least one vertex
   * @param[in] eleMap  : 'Map' from the grid-global cell ID to its index within 'eles' vector (or -1 if blanked cell)
   * @param[in] centroid: Centriod of current grid partition
   * @param[in] extents : x,y,z extents (max-min) of current grid partition
   */
  void matchOversetPoints(vector<shared_ptr<ele>> &eles, const vector<int> &eleMap, const point &minPt = point({0,0,0}), const point &maxPt = point({0,0,0}));

  /*!
   * \brief Setup all communication for unblanked cells and faces
   *
   * For each unblanked cell, setup a superMesh from donor elements on this grid
   * --The quadrature order used on the superMesh will be order quadOrder
   * For each blanked face, remove its points from the communicator
   * For each unblanked face, add its points to the communicator
   *
   */
  void matchUnblankCells(vector<shared_ptr<ele>> &eles, unordered_set<int>& unblankCells, vector<int>& eleMap, int quadOrder);

  void performGalerkinProjection(vector<shared_ptr<ele> >& eles, map<int, oper>& opers, vector<int>& eleMap, int order);

  void performProjection_static(vector<shared_ptr<ele> >& eles, vector<int>& eleMap, int order);

  //! Integrate the solution error over the entire domain, accounting for overset overlap
  vector<double> integrateErrOverset(vector<shared_ptr<ele> >& eles, map<int, oper>& opers, vector<int>& iblankCell, vector<int>& eleMap, int order, int quadOrder);

  //! Perform the interpolation and communicate data across all grids
  void exchangeOversetData(vector<shared_ptr<ele>> &eles, map<int, oper>& opers, vector<int> &eleMap);

  //! Perform the interpolation and communicate gradient across all grids
  void exchangeOversetGradient(vector<shared_ptr<ele>> &eles, map<int, oper>& opers, vector<int> &eleMap);

  /*!
   * \brief Gather a distributed dataset so that every rank has the full, organized dataset
   *
   * @param[in] nPieces      : Number of pieces of data contributed by current rank
   * @param[in] stride       : Number of [typename T] values for each piece of data
   * @param[in] values       : Values being contributed from current rank
   * @param[out] nPieces_rank: Number of pieces of data for each rank
   * @param[out] values_all  : The collected and organized values for all grids/ranks
   */
  template<typename T>
  void gatherData(int nPieces, int stride, T *values, vector<int>& nPieces_rank, vector<T> &values_all);

  /*!
   * \brief Distribute a scattered dataset to the proper grid/rank
   *
   * @param[in] nPiecesSend : Number of pieces of data this rank is sending to each grid
   * @param[in] nPiecesRecv : Number of pieces of data this rank will be receiving from each grid
   * @param[in] sendInds    : Destination indices on receiving rank for data being sent (range 0:nPieces_rank-1)
   * @param[in/out] recvInds: Destination indices on current grid for data being received (range 0:nPieces-1)
   * @param[in] nPieces_rank: Number of pieces of data for each rank on current grid
   * @param[in] values_send : Dataset being sent to each grid (size nPiecesSend x stride)
   * @param[out] values_recv: Dataset received from all other grids
   * @param[in] stride      : Number of [typename T] values for each piece of data
   * @param[in] matchInds   : Whether each rank needs to receive recvInds for its data or not
   */
  template<typename T>
  void sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<vector<int>> &recvInds, vector<matrix<T>> &sendVals, matrix<T> &recvVals, int stride, bool matchInds);

  template<typename T>
  void sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<vector<int>> &recvInds, vector<matrix<T>> &sendVals, vector<matrix<T>> &recvVals, int stride);

  template<typename T>
  void sendRecvData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<matrix<T>> &sendVals, vector<matrix<T>> &recvVals, int stride);

  //! Using nPiecesIn, resize nPiecesOut for each rank
  void setupNPieces(vector<int> &nPiecesIn, vector<int> &nPiecesOut);

  //! Setup the list of overset-boundary flux points to interpolate data to
  void setupOverFacePoints(vector<shared_ptr<overFace> >& overFaces);

  //! Setup the list of fringe/receptor solution points to interpolate data to
  void setupFringeCellPoints(vector<shared_ptr<ele> >& eles, const unordered_set<int>& fringeCells, const vector<int>& eleMap);

  //! Transfer the exchanged U_in data to the fringe cells
  void transferEleData(vector<shared_ptr<ele> >& eles, const unordered_set<int>& fringeCells, const vector<int>& eleMap);

private:

  /* ---- For Static Cases using Field Interpolation ---- */
  vector<matrix<double>> qpts, qptsD_ref, donorBasis, massMatTDRow, ubLHS;
  vector<vector<int>> targetID, donorID;
  vector<vector<int>> recvInds;

  //! For use with ADT in 2D
  vector<int> eleList;
};
