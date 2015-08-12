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
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2015 Jacob Crabill.
 *
 */
#pragma once

#include "global.hpp"

#include <map>
#include <set>
#include <vector>

class oper;
class overFace;

#include "ele.hpp"
#include "input.hpp"
#include "operators.hpp"
#include "overFace.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "mpi.h"
#endif

struct dataExchange
{
  vector<int> nPts_rank;           //! Number of fringe points for each rank of current grid
  vector<vector<int>> foundPts;    //! IDs of receptor points from each grid which were found to lie within current grid
  vector<vector<int>> foundRank;   //! gridRank of this process for each found point (for benefit of other processes; probably not needed)
  vector<vector<int>> foundEles;   //! Ele ID which each matched point was found to lie within
  vector<vector<point>> foundLocs; //! Reference location within donor ele of each matched receptor point

  int nOverPts;                 //! Number of overset (receptor) points on this grid
  matrix<double> overPts;       //! Physical locations of the receptor points on this grid
  vector<int> nPtsRecv;         //! Number of points incoming from each grid (across interComm)
  vector<int> nPtsSend;         //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvPts;  //! Point IDs which will be received from each grid (across interComm) (counter to foundPts)

  matrix<double> U_in;          //! Data received from other grid(s)
  vector<matrix<double>> U_out; //! Interpolated data being sent to other grid(s)
};

class overComm
{
public:
  overComm();

  int nGrids;
  int nprocPerGrid;
  int gridID;
  int gridRank;

  int nFields;

  input *params;  //! Simulation input parameters

#ifndef _NO_MPI
  MPI_Comm gridComm;
  MPI_Comm interComm;
#endif

  /* --- Variables for Exchanging Data at Overset Faces --- */

  vector<int> nPts_rank;           //! Number of fringe points for each rank of current grid
  vector<vector<int>> foundPts;    //! IDs of receptor points from each grid which were found to lie within current grid
  vector<vector<int>> foundRank;   //! gridRank of this process for each found point (for benefit of other processes; probably not needed)
  vector<vector<int>> foundEles;   //! Ele ID which each matched point was found to lie within
  vector<vector<point>> foundLocs; //! Reference location within donor ele of each matched receptor point

  int nOverPts;                 //! Number of overset (receptor) points on this grid
  matrix<double> overPts;       //! Physical locations of the receptor points on this grid
  vector<int> nPtsRecv;         //! Number of points incoming from each grid (across interComm)
  vector<int> nPtsSend;         //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvPts;  //! Point IDs which will be received from each grid (across interComm) (counter to foundPts)

  matrix<double> U_in;          //! Data received from other grid(s)
  vector<matrix<double>> U_out; //! Interpolated data being sent to other grid(s)

  /* --- Variables for Exchanging Data on Unblanked Cells --- */

  vector<int> nCells_rank;              //! Number of unblanked cells for each rank of current grid
  vector<vector<int>> foundCells;       //! IDs of unblanked cells from each grid which were found to overlap this grid
  vector<matrix<int>> foundCellDonors;  //! Donor cells on this grid for each found unblanked cell
  vector<vector<int>> foundCellNDonors; //! Number of donor cells on this grid for each found unblanked cell

  int nUnblanks;          //! Number of unblank cells on this grid
  vector<int> unblanks;   //! Cells from this grid which need to be unblanked
  vector<int> nCellsRecv;         //! Number of points incoming from each grid (across interComm)
  vector<int> nCellsSend;         //! Number of points outgoing to each grid (across interComm)
  vector<vector<int>> recvCells;  //! Cell IDs which will be received from each grid (across interComm) (counter to foundPts)

  //! Local supermesh of donor elements for each cell needing to be unblanked
  vector<superMesh> donors;

  /* --- Member Functions --- */

  void setup(input *_params, int _nGrids, int _gridID, int _gridRank, int _nprocPerGrid);

  //! Match up each overset-face flux point to its donor grid and element
  void matchOversetPoints(vector<ele>& eles, vector<shared_ptr<overFace> >& overFaces);

  /*!
   * \brief Setup all communication for unblanked cells and faces
   *
   * For each unblanked cell, setup a superMesh from donor elements on this grid
   * --The quadrature order used on the superMesh will be order quadOrder
   * For each blanked face, remove its points from the communicator
   * For each unblanked face, add its points to the communicator
   *
   */
  void matchOversetUnblanks(vector<ele> &eles, vector<shared_ptr<overFace> >& overFaces, set<int>& unblankCells,
                            set<int> &blankedOFaces, set<int>& unblankOFaces, vector<int>& eleMap, vector<int>& faceMap, int quadOrder);

  //! Perform the interpolation and communicate data across all grids
  void exchangeOversetData(vector<ele> &eles, map<int, map<int,oper> > &opers);

  /*!
   * \brief Gather a distributed dataset so that every rank has the full, organized dataset
   *
   * @param[in] nPieces      : Number of pieces of data contributed by current rank
   * @param[in] stride       : Number of [typename T] values for each piece of data
   * @param[in] values       : Values being contributed from current rank
   * @param[out] nPieces_rank: Number of pieces of data for each rank on current grid
   * @param[out] nPieces_grid: Number of pieces of data for each grid
   * @param[out] values_alll : The collected and organized values for all grids/ranks
   */
  template<typename T>
  void gatherData(int nPieces, int stride, T *values, vector<int> &nPieces_rank, vector<int> &nPieces_grid, vector<T> &values_all);

  /*!
   * \brief Distribute a scattered dataset to the proper grid/rank
   *
   * @param[in] nPiecesSend : Number of pieces of data this rank is sending to each grid
   * @param[in] nPiecesRecv : Number of pieces of data this rank will be receiving from each grid
   * @param[in] sendInds    : Destination indices on receiving grid for data being sent (range 0:nPieces_grid-1)
   * @param[in] nPieces_rank: Number of pieces of data for each rank on current grid
   * @param[in] values_send : Dataset being sent to each grid (size nPiecesSend x stride)
   * @param[out] values_recv: Dataset received from all other grids
   * @param[in] stride      : Number of [typename T] values for each piece of data
   */
  template<typename T>
  void distributeData(vector<int> &nPiecesSend, vector<int> &nPiecesRecv, vector<vector<int>> &sendInds, vector<int> &nPieces_rank, vector<matrix<T>> &values_send, matrix<T> &values_recv, int stride);

private:
#ifndef _NO_MPI
  template<typename T> MPI_Datatype getMpiDatatype(void);
#endif
};
