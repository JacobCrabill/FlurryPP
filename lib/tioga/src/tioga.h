/**
 * Topology Indpendent Overset Grid Assembler (TIOGA)
 * Base class and dependencies
 * The methods of this class are invoked from tiogaInterface.C
 *
 *  Jay Sitaraman 02/24/2014
 */

#pragma once

#include <set>

class MeshBlock;

#include "parallelComm.h"
#include "solver.hpp"

class tioga
{
 private :
  int nblocks;
  MeshBlock *mb;
  int nmesh;
  HOLEMAP *holeMap;
  MPI_Comm scomm;
  parallelComm *pc;
  int isym;
  int ierr;
  int mytag;
  int myid,numprocs;
  int *sendCount;
  int *recvCount;
  OBB *obblist;
  int iorphanPrint;

  solver *Solver;

 public:
  int ihigh;

  /** basic constuctor */
  tioga();

  /** basic destructor */
  ~tioga();

  /** set communicator */
  void setCommunicator(MPI_Comm communicator,int id_proc,int nprocs);

  /** registerGrid data */

  void registerGridData(int btag,int nnodes,double *xyz,int *ibl, int nwbc,int nobc,
			       int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn);

  void profile(void);

  int findPointDonor(double* x_pt);

  std::unordered_set<int> findCellDonors(double* bbox);

  void exchangeBoxes(void);

  void exchangeSearchData(void);

  void exchangePointSearchData(void);

  void exchangeDonors(void);

  /** perform overset grid connectivity */

  void performConnectivity(void);

  void performConnectivityHighOrder(void);

  /** update data */

  void dataUpdate(int nvar,double *q,int interptype) ;

  void dataUpdate_highorder(int nvar,double *q,int interptype) ;

  /** get hole map for each mesh */

  void getHoleMap(void);

  /** output HoleMaps */

  void outputHoleMap(void);

  void writeData(int nvar,double *q,int interptype);

  void getDonorCount(int *dcount, int *fcount);

  void getDonorInfo(int *receptors,int *indices,double *frac,int *dcount);

  /** set symmetry bc */
  void setSymmetry(int syminput) { isym=syminput;}

  /** set resolutions for nodes and cells */
  void setResolutions(double *nres,double *cres);

  void set_cell_iblank(int *iblank_cell);

  //! Set callback functions for high-order processing. See MeshBlock.h for details
  void setcallback(void (*f1)(int*, int*),
		    void (*f2)(int *,int *,double *),
		    void (*f3)(int *,double *,int *,double *),
		    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
        void (*f5)(int *,int *,double *,int *,int*,double *));

  void setcallback(solver* _solver);
};




