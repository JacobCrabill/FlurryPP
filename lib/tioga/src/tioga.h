/**
 * Topology Indpendent Overset Grid Assembler (TIOGA)
 * Base class and dependencies
 * The methods of this class are invoked from tiogaInterface.C
 *
 *  Jay Sitaraman 02/24/2014 
 */
#include "MeshBlock.h"
#include "parallelComm.h"

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

 public:
  int ihigh;

  /** basic constuctor */
  tioga() {
    mb = NULL;
    holeMap=NULL;
    pc=NULL;
    sendCount=NULL;
    recvCount=NULL;
    obblist=NULL;
    isym=2;
    ihigh=0;
  };
 
  /** basic destructor */
  ~tioga(); 
  
  /** set communicator */
  void setCommunicator(MPI_Comm communicator,int id_proc,int nprocs);

  /** registerGrid data */

  void registerGridData(int btag,int nnodes,double *xyz,int *ibl, int nwbc,int nobc,
			       int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn);

  void profile(void);

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
  void setResolutions(double *nres,double *cres) { mb->setResolutions(nres,cres);}

  void set_cell_iblank(int *iblank_cell)
  {
   mb->set_cell_iblank(iblank_cell);
  }

  void setcallback(void (*f1)(int*, int*),
		    void (*f2)(int *,int *,double *),
		    void (*f3)(int *,double *,int *,double *),
		    void (*f4)(int *,double *,int *,int *,double *,double *,int *),
		   void (*f5)(int *,int *,double *,int *,int*,double *))
  {
    mb->setcallback(f1,f2,f3,f4,f5);
    ihigh=1;
  }
};
      
  


