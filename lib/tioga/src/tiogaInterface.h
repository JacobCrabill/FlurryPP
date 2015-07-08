#ifndef TIOGAINTERFACE_H
#define TIOGAINTERFACE_H

//
// All the interfaces that are
// accessible to third party f90 and C
// flow solvers
//
//
// Jay Sitaraman
// 02/24/2014
//
// Note that Fortran requires the trailing underscore in
// C function names, but will be called without it
// i.e. in Fortran, use: CALL TIOGA_INIT(...)
extern "C" {
  void tioga_init_(int* scomm);

  void tioga_registergrid_data_(int* btag,int* nnodes,double* xyz,int* ibl,int* nwbc, int* nobc,int* wbcnode,
                                int* obcnode,int* ntypes,...);

  void tioga_preprocess_grids_(void);

  void tioga_performconnectivity_(void);

  void tioga_performconnectivity_highorder_(void);

  void tioga_dataupdate_(double* q,int* nvar,char* itype);

  void tioga_writeoutputfiles_(double* q,int* nvar,char* itype);

  void tioga_getdonorcount_(int* dcount,int* fcount);

  void tioga_getdonorinfo_(int* receptors,int* indices,double* frac,int* dcount);

  void tioga_setsymmetry_(int* isym);

  void tioga_setresolutions_(double* nres,double* cres);

  void tioga_setcelliblank_(int* iblank_cell);

  /*!
   * Assign callback functions for high-order method processing
   * f1: get_nodes_per_cell (cellID, nodesPerCell)
   * f2: get_recepter_nodes (cellID, pointsPerCell, receptor_nodes_XYZ)
   * f3: donor_inclusion_test (cellID, posXYZ, passFlag, rst[ipoint])
   * f4: donor_frac (cellID, xsearch, nweights, inode, frac, rst, ndim)
   * f5: convert_to_modal (cellID, pointsPerCell, nodalVals, nPts, indexOut, modalVals)
   */
  void tioga_set_highorder_callback_(void (*f1)(int*, int*),
                                     void (*f2)(int *,int *,double *),
                                     void (*f3)(int *,double *,int *,double *),
                                     void (*f4)(int *,double *,int *,int *,double *,double *,int *),
                                     void (*f5)(int *,int *,double *,int *,int *,double *));


  void tioga_delete_(void);
}

#endif // TIOGAINTERFACE_H
