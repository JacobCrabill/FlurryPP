/**
 * MeshBlock class - container and functions for generic unstructured grid partition in 3D
 *         
 * Jay Sitaraman
 * 02/20/2014
 */
#include "codetypes.h"
#include "ADT.h"

class MeshBlock
{
 private:
  int nnodes;  /** < number of grid nodes */
  int ncells;  /** < total number of cells */
  int ntypes;  /** < number of different types of cells */
  int *nv;     /** < number of vertices for each types of cell */
  int *nc;     /**  < number of each of different kinds of cells (tets, prism, pyramids, hex etc) */
  int nobc;    /** < number of overset boundary nodes */
  int nwbc;    /** < number of wall boundary nodes */

  double *x;        /** < grid nodes x[3*nnodes] */
  int *iblank;      /** < iblank value for each grid node */
  int *iblank_cell; /** < iblank value at each grid cell */

  int **vconn;        /** < connectivity of each kind of cell */
  int *wbcnode;     /** < wall boundary node indices */
  int *obcnode;     /** < overset boundary node indices */

  double *nodeRes;  /** < node resolution  */
  double *userSpecifiedNodeRes;
  double *userSpecifiedCellRes;
  double *elementBbox; /** < bounding box of the elements */
  int *elementList;    /** < list of elements in */

  //! Alternating digital tree library
  ADT *adt;   /** < Digital tree for searching this block */

  DONORLIST **donorList;      /**< list of donors for the nodes of this mesh */

  int ninterp;              /**< number of interpolations to be performed */
  INTERPLIST *interpList;   /**< list of donor nodes in my grid, with fractions and information of
                                 who they donate to */ 
  int *interp2donor;

  INTEGERLIST *cancelList;  /** receptors that need to be cancelled because of */
  int ncancel;              /** conflicts with the state of their donors */

  /* ---- Callback functions for high-order overset connectivity ---- */

  /*!
   * \brief Get the number of solution points in given cell
   */
  void (*get_nodes_per_cell)(int* cellID, int* nNodes);

  /*!
   * \brief Get the physical position of solution points in given cell
   *
   * input: cellID, nNodes
   * output: xyz [size: nNodes x 3, row-major]
   */
  void (*get_receptor_nodes)(int* cellID, int* nNodes, double* xyz);

  /*!
   * \brief Determine whether a point (x,y,z) lies within a cell
   *
   * Given a point's physical position, determine if it is contained in the
   * given cell; if so, return the reference coordinates of the point
   *
   * @param[in]  cellID    ID of cell within current mesh
   * @param[in]  xyz       Physical position of point to test
   * @param[out] passFlag  Is the point inside the cell? (no:0, yes:1)
   * @param[out] rst       Position of point within cell in reference coordinates
   */
  void (*donor_inclusion_test)(int* cellID, double* xyz, int* passFlag, double* rst);

  /*!
   * \brief Get interpolation points & weights for current cell,
   *        given a point in reference coordinates
   *
   * @param[in]  cellID    ID of cell within current mesh
   * @param[in]  xyz       Physical position of receptor point
   * @param[out] nweights  Number of interpolation points/weights to be used
   * @param[out] inode     Indices of donor points within global solution array
   * @param[out] weights   Interpolation weights for each donor point
   * @param[in]  rst       Reference coordinates of receptor point within cell
   * @param[in]  ndim      Amount of memory allocated to 'frac' (# of doubles)
   */
  void (*donor_frac)(int* cellID, double* xyz, int* nweights, int* inode, double* weights, double* rst, int* ndim);

  void (*convert_to_modal)(int *,int *,double *,int *,int *,double *);

  /* ---- End High-Order Callback Function Definitions ---- */

  int nreceptorCells;      /** number of receptor cells */
  int *ctag;               /** index of receptor cells */
  int *pointsPerCell;      /** number of receptor points per cell */
  int maxPointsPerCell;    /** max of pointsPerCell vector */
  double *rxyz;            /**  point coordinates */
  int ipoint; 

 public :
  int nfringe;
  int meshtag; /** < tag of the mesh that this block belongs to */

  //! oriented bounding box of this partition
  OBB *obb;

  int nsearch;        /** < number of query points to search in this block */
  int *isearch;       /** < index of query points in the remote process */
  double *xsearch;    /** < coordinates of the query points */
  double *rst;        /** < natrural coordinates */
  int *donorId;       /** < donor (cell?) indices for those found */
  int donorCount;
  int myid;
  double *cellRes;    /** < resolution for each cell */
  int ntotalPoints;   /**  total number of extra points to interpolate */
  int ihigh;
  int ninterp2;            /** < number of interpolants for high-order points */
  INTERPLIST *interpList2; /** < list for high-interpolation points */
  
  /** basic constructor */
  MeshBlock() {
    nv=NULL;
    nc=NULL;
    x=NULL;
    iblank=NULL;
    iblank_cell=NULL;
    vconn=NULL;
    wbcnode=NULL;
    obcnode=NULL;
    cellRes=NULL;
    nodeRes=NULL;
    elementBbox=NULL;
    elementList=NULL;
    adt=NULL;
    donorList=NULL;
    interpList=NULL;
    interp2donor=NULL;
    obb=NULL;
    nsearch=0;
    isearch=NULL;
    xsearch=NULL;
    donorId=NULL;
    adt=NULL;
    cancelList=NULL;
    userSpecifiedNodeRes=NULL;
    userSpecifiedCellRes=NULL;
    nfringe=2;
    // new vars
    ctag=NULL;
    pointsPerCell=NULL;
    maxPointsPerCell=0;
    rxyz=NULL;
    ntotalPoints=0;
    rst=NULL;
    ihigh=0;
    ipoint=0;
    interpList2=NULL;
  }

  /** basic destructor */
  ~MeshBlock();
      
  void preprocess(void);

  void tagBoundary(void);
  
  void writeGridFile(int bid);

  void writeFlowFile(int bid,double *q,int nvar,int type);
  
  void setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci, 
	       int *wbcnodei,int *obcnodei,
	       int ntypesi, int *nvi, int *nci, int **vconni);

  void setResolutions(double *nres,double *cres);    
	       
  void search();

  void writeOBB(int bid);

  void updateSolnData(int inode,double *qvar,double *q,int nvar,int interptype);

  int getNinterp(void) {return ninterp;}

  void getInterpolatedSolution(int *nints,int *nreals,int **intData,double **realData,double *q,
			       int nvar, int interptype);
  
  void checkContainment(int *cellIndex,int adtElement,double *xsearch);

  void getWallBounds(int *mtag,int *existWall, double wbox[6]);
  
  void markWallBoundary(int *sam,int nx[3],double extents[6]);

  void getQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,
		      double **realData);
  

  /** routines that do book keeping */

  void getDonorPacket(PACKET *sndPack, int nsend);

  void initializeDonorList();
  
  void insertAndSort(int pointid,int senderid,int meshtag, int remoteid, double donorRes);
  
  void processDonors(HOLEMAP *holemap, int nmesh,int **donorRecords,double **receptorResolution,
		     int *nrecords);

  void initializeInterpList(int ninterp_input);
  
  void findInterpData(int* recid,int irecord, double receptorRes);

  void set_ninterp(int);

  void getCancellationData(int *nints, int **intData);

  void cancelDonor(int irecord);

  void getInterpData(int *nrecords, int **intData);

  void clearIblanks(void);

  void setIblanks(int inode);

  void getDonorCount(int *dcount,int *fcount);

  void getDonorInfo(int *receptors,int *indices, double *frac);
  //
  // routines for high order connectivity and interpolation
  //
  void getCellIblanks(void);
  void set_cell_iblank(int *iblank_cell_input)
  {
    iblank_cell=iblank_cell_input;
  }

  void setcallback(void (*f1)(int*, int*),
                   void (*f2)(int *,int *,double *),
                   void (*f3)(int *,double *,int *,double *),
                   void (*f4)(int *,double *,int *,int *,double *,double *,int *),
                   void (*f5)(int *,int *,double *,int *,int*,double *))
  {
    get_nodes_per_cell=f1;    // (cellID, pointsPerCell*)                       [* = output]
    get_receptor_nodes=f2;    // (cellID, pointsPerCell, receptor_nodes_XYZ*)
    donor_inclusion_test=f3;  // (cellID, posXYZ, passFlag*, rst[iPoint]*)
    donor_frac=f4;            // (cellID, xsearch, nweights, inode, frac*, rst, ndim)
    convert_to_modal=f5;      // (cellID, pointsPerCell, nodalVals, nPts, indexOut*, modalVals*)
  }

  void writeCellFile(int);

  void getInternalNodes(void);

  void getExtraQueryPoints(OBB *obb,int *nints,int **intData,int *nreals,double **realData);

  void processPointDonors(void);

  void getInterpolatedSolutionAtPoints(int *nints,int *nreals,int **intData,double **realData,
                                       double *q,int nvar, int interptype);

  void updatePointData(double *q,double *qtmp,int nvar,int interptype);

  void outputOrphan(FILE *fp,int i)
  {
    fprintf(fp,"%f %f %f\n",rxyz[3*i],rxyz[3*i+1],rxyz[3*i+2]);
  }
};
