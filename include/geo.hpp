/*!
 * \file geo.hpp
 * \brief Header file for geometry class
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
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "global.hpp"

#ifndef _NO_MPI
class tioga;
#endif

#include "ele.hpp"
#include "input.hpp"
#include "face.hpp"
#include "mpiFace.hpp"
#include "overFace.hpp"
#include "solver.hpp"
#include "superMesh.hpp"

#ifndef _NO_MPI
#include "tioga.h"
#endif

#define NORMAL  1
#define HOLE    0
#define FRINGE -1

class geo
{
public:
  geo();

  ~geo();

  /* === Primay setup routines === */

  //! Setup the geomery using input parameters
  void setup(input* params);

  //! Take the basic connectivity data and generate the rest
  void processConnectivity();

  //! Create the elements and faces needed for the simulation
  void setupElesFaces(vector<ele> &eles, vector<shared_ptr<face> > &faces, vector<shared_ptr<mpiFace> > &mpiFacesVec, vector<shared_ptr<overFace> >& overFacesVec);

  //! Update nodal positions and velocities for moving-grid cases
  void moveMesh(void);

  /* === Helper Routines === */

  //! Read essential connectivity from a Gmsh mesh file
  void readGmsh(string fileName);

  //! Create a simple Cartesian mesh from input parameters
  void createMesh();

  //! Update connectivity / node-blanking for overset grids
  void registerGridDataTIOGA();

  /*!
   * \brief Call TIOGA to re-process overset connectivity
   *
   * Called once during pre-processing by default; re-call each iteration
   * for moving-mesh cases
   */
  void updateOversetConnectivity();

  //! Have TIOGA output the mesh along with nodal IBLANK values
  void writeOversetConnectivity();

  /* ---- My Overset Functions ---- */

  //! Setup the 'connectivity' between the overset interpolation points & donor grids/cells
  void matchOversetPoints(vector<ele> &eles, struct dataExchange& exchange);

  //! Send / Receive interpolated data to proper grid and rank
  void exchangeOversetData(struct dataExchange &exchange);

  void matchOversetDonors(vector<ele> &eles, vector<superMesh> &donors);

  void setupUnblankElesFaces(vector<ele> &eles, vector<shared_ptr<face>> &faces, vector<shared_ptr<mpiFace>> &mpiFacesVec, vector<shared_ptr<overFace>> &overFacesVec);

  int nDims, nFields;
  int nEles, nVerts, nEdges, nFaces, nIntFaces, nBndFaces, nMpiFaces, nOverFaces;
  int nBounds;  //! Number of boundaries
  int meshType;

  // Basic [essential] Connectivity Data
  matrix<int> c2v;
  matrix<double> xv;      //! Current physical position of vertices [static or moving grids]
  vector<point> xv_new;   //! Physical position of vertices for next time step [moving grids]
  vector<point> xv0;      //! Initial position of vertices [moving grids]
  matrix<double> gridVel; //! Grid velocity of vertices

  // Additional Connectivity Data
  matrix<int> c2e, c2b, e2c, e2v, v2e, v2v, v2c;
  matrix<int> c2f, f2v, f2c;
  vector<int> v2nv, v2nc, c2nv, c2nf, f2nv, ctype;
  vector<int> intFaces, bndFaces, mpiFaces, overFaces, mpiCells;
  vector<int> bcList;            //! List of boundary conditions for each boundary
  vector<int> bcType;            //! Boundary condition for each boundary edge
  matrix<int> bndPts;            //! List of node IDs on each boundary
  vector<int> nBndPts;           //! Number of points on each boudary
  vector<matrix<int> > bcFaces;  //! List of nodes on each face (edge) for each boundary condition
  vector<int> nFacesPerBnd;      //! List of # of faces on each boundary
  vector<int> procR;             //! What processor lies to the 'right' of this face
  vector<int> faceID_R;            //! The local mpiFace ID of each mpiFace on the opposite processor
  vector<int> gIC_R;             //! The global cell ID of the right cell on the opposite processor
  vector<int> mpiLocF;           //! Element-local face ID of MPI Face in left cell
  vector<int> mpiLocF_R;         //! Element-local face ID of MPI Face in right cell
  vector<int> mpiPeriodic;       //! Flag for whether an MPI face is also a periodic face
  vector<bool> isBnd; // might want to change this to "int" and have it store WHICH boundary the face is on (-1 for internal)

  /* --- Overset-Related Variables --- */
  int nGrids;             //! Number of distinct overset grids
  int nprocPerGrid;       //! Number of MPI processes assigned to each (overset) grid block
  int gridID;             //! Which (overset) grid block is this process handling
  int gridRank;           //! MPI rank of process *within* the grid block [0 to nprocPerGrid-1]
  vector<int> iblank;     //! Output of TIOGA: flag for whether vertex is normal, blanked, or receptor
  vector<int> iblankCell; //! Output? of TIOGA: flag for whether cell is normal, blanked, or receptor
  vector<int> iblankFace; //! Flag for whether a face is normal, blanked, or receptor
  vector<int> iwall;      //! List of nodes on wall boundaries
  vector<int> iover;      //! List of nodes on overset boundaries

  vector<int> eleMap;     //! For overset meshes where some cells are blanked, map from 'ic' to 'eles' index
  vector<int> faceMap;
  //vector<int> bfaceMap;
  //vector<int> mFaceMap;
  //vector<int> oFaceMap;

#ifndef _NO_MPI
  MPI_Comm gridComm;  //! Intra-grid communicator
  MPI_Comm interComm; //! Inter-grid communicator (matched by gridRank)
#endif

  /* --- Moving-Overset-Grid-Related Variables --- */
  set<int> holeCells;     //! List of cells in mesh which are currently blanked
  set<int> holeFaces;     //! List of cells in mesh which are currently blanked
  set<int> unblankCells;  //! List of non-existing cells which, due to motion, must be un-blanked
  set<int> unblankFaces;  //! List of non-existing faces which, due to motion, must be un-blanked
  set<int> unblankOFaces; //! List of non-existing faces which, due to motion, must be un-blanked as overset faces
  set<int> blankCells;    //! List of existing cells which, due to motion, must be blanked
  set<int> blankFaces;    //! List of existing faces which, due to motion, must be blanked

  // Outgoing (interpolated) data
  vector<int> interpProc;       //! Processor ID for all interpolation points
  vector<int> interpCell;       //! For each interpolation point, cell which it lies within
  matrix<double> interpPtsPhys; //! Physical positions of fringe points on other grids to interpolate to
  matrix<double> interpPtsRef;  //! Reference position (within interpCell) of fringe points on other grids to interpolate to

  // Incoming (overset) data
//  vector<int> overProc;         //! Donor-processor ID for each fringe point
//  matrix<double> overPtsPhys;   //! Physical positions of each fringe point
//  vector<point> overPts;        //! Physical positions of fringe points

#ifndef _NO_MPI
  shared_ptr<tioga> tg;  //! Pointer to Tioga object for processing overset grids
#endif
  int* nodesPerCell;   //! Pointer for Tioga to know # of nodes for each element type
  array<int*,1> conn;  //! Pointer to c2v for each element type [but only 1, so will be size(1)]
  matrix<int> tg_c2v;  //! 'Cleaned' c2v for Tioga (when quadratic elements present, normal c2v won't work)

  // Data for communication between grids
//  vector<int> nPts_rank;           //! Number of fringe points for each rank of current grid
//  vector<vector<int>> foundPts;    //! IDs of receptor points from each grid which were found to lie within current grid
//  vector<vector<int>> foundRank;   //! gridRank of this process for each found point (for benefit of other processes; probably not needed)
//  vector<vector<int>> foundEles;   //! Ele ID which each matched point was found to lie within
//  vector<vector<point>> foundLocs; //! Reference location within ele of each matched receptor point

//  vector<int> nPtsRecv; //! Number of points incoming from each grid (across interComm)
//  vector<int> nPtsSend; //! Number of points outgoing to each grid (across interComm)
//  vector<vector<int>> recvPts; //! Point IDs which will be received from each grid (across interComm) (counter to foundPts)

private:

  input *params;

  /* --- MPI-Related Varialbes (global vs. local data) --- */
  matrix<int> c2v_g;    //! Global element connectivity
  matrix<double> xv_g;  //! Global mesh node locations
  vector<int> ic2icg;   //! Local cell to global cell index
  vector<int> iv2ivg;   //! Local vertex to global vertex index
  vector<int> ctype_g, c2ne_g, c2nv_g; //! Global element info
  matrix<int> bndPts_g;  //! Global lists of points on boundaries
  vector<int> nBndPts_g; //! Global number of points on each boundary
  map<int,int> bcIdMap;  //! Map from Gmsh boundary ID to Flurry BC ID
  int nEles_g, nVerts_g;

  void processConn2D(void);
  void processConn3D(void);

  //! Using Tioga's nodal iblanks, set iblank values for all cells and faces
  void setCellFaceIblanks();

  //! Match up pairs of periodic boundary faces
  void processPeriodicBoundaries(void);

  //! Check if two given periodic edges match up
  bool checkPeriodicFaces(int *edge1, int *edge2);
  bool checkPeriodicFaces3D(vector<int> &face1, vector<int> &face2);
  bool comparePeriodicMPI(vector<int> &face1, vector<int> &face2);

  //! Compare the orientation (rotation in ref. space) betwen the local faces of 2 elements
  int compareOrientation(int ic1, int ic2, int f1, int f2);

  //! Compare the orientation (rotation in ref. space) betwen the local faces of 2 elements across MPI boundary
  int compareOrientationMPI(int ic1, int ic2, int f1, int f2, int isPeriodic);

  //! For MPI runs, partition the mesh across all processors
  void partitionMesh(void);

  //! For MPI runs, match internal faces across MPI boundaries
  void matchMPIFaces();

  //! Compare two faces [lists of nodes] to see if they match [used for MPI]
  bool compareFaces(vector<int> &face1, vector<int> &face2);

};
