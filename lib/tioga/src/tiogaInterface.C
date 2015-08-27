#include "tiogaInterface.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "tioga.h"
#include "globals.h"

#include "solver.hpp" // Include Flurry's 'solver' class for callback functions

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
  void tioga_init_(int* scomm)
  {
    int id_proc,nprocs;
    MPI_Comm tcomm;
    //tcomm=(MPI_Comm) (*scomm);
    tcomm=MPI_COMM_WORLD;
    //
    tg=new tioga[1];
    //
    //MPI_Comm_rank(MPI_COMM_WORLD,&id_proc);
    //MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(tcomm,&id_proc);
    MPI_Comm_size(tcomm,&nprocs);
    //
    tg->setCommunicator(tcomm,id_proc,nprocs);
    nc=NULL;
    nv=NULL;
    vconn=NULL;
  }


  void tioga_registergrid_data_(int* btag,int* nnodes,double* xyz,int* ibl,int* nwbc, int* nobc,int* wbcnode,
                                int* obcnode,int* ntypes,...)
  {
    va_list arguments;
    int i;

    va_start(arguments, *ntypes);
    nv=(int *) malloc(sizeof(int)*(*ntypes));
    nc=(int *) malloc(sizeof(int)*(*ntypes));
    vconn=(int **)malloc(sizeof(int *)*(*ntypes));
    for(i=0; i<*ntypes; i++)
    {
      nv[i]=*(va_arg(arguments, int *));
      nc[i]=*(va_arg(arguments, int *));
      vconn[i]=va_arg(arguments, int *);
    }
    tg->registerGridData(*btag,*nnodes,xyz,ibl,*nwbc,*nobc,wbcnode,obcnode,*ntypes,nv,nc,vconn);
  }

  void tioga_preprocess_grids_(void)
  {
    tg->profile();
  }

  void tioga_performconnectivity_(void)
  {
    tg->performConnectivity();
  }

  void tioga_performconnectivity_highorder_(void)
  {
    tg->performConnectivityHighOrder();
  }

  void tioga_dataupdate_(double* q,int* nvar,int* itype)
  {
    // 0: row, 1: column
    interptype = itype;

    if (itype != 0 && itype != 1)
    {
      printf("#tiogaInterface.C:dataupdate_:unknown data orientation\n");
      return;
    }

    if (tg->ihigh==0)
    {
      tg->dataUpdate(*nvar,q,interptype);
    }
    else
    {
      tg->dataUpdate_highorder(*nvar,q,interptype);
    }
  }

  void tioga_writeoutputfiles_old_(double* q,int* nvar,char* itype)
  {
    char itypeC[4] = {itype[0],itype[1],itype[2],'\0'};
    int interptype;
    if (strcmp(itypeC,"row")==0)
    {
      interptype=0;
    }
    else if (strcmp(itypeC,"col")==0)
    {
      interptype=1;
    }
    else
    {
      printf("#tiogaInterface.C:dataupdate_:unknown data orientation %s\n",itype);
      return;
    }
    tg->writeData(*nvar,q,interptype);
  }

  void tioga_writeoutputfiles_(double* q,int* nvar,int* itype)
  {
    int interptype = *itype;

    if (interptype !=0 && interptype != 1)
    {
      printf("#tiogaInterface.C:dataupdate_:unknown data orientation %s\n",itype);
      return;
    }
    tg->writeData(*nvar,q,interptype);
  }

  void tioga_getdonorcount_(int* dcount,int* fcount)
  {
    tg->getDonorCount(dcount,fcount);
  }
  void tioga_getdonorinfo_(int* receptors,int* indices,double* frac,int* dcount)
  {
    tg->getDonorInfo(receptors,indices,frac,dcount);
  }

  void tioga_setsymmetry_(int* isym)
  {
    tg->setSymmetry(*isym);
  }

  void tioga_setresolutions_(double* nres,double* cres)
  {
    tg->setResolutions(nres,cres);
  }

  void tioga_setcelliblank_(int* iblank_cell)
  {
    tg->set_cell_iblank(iblank_cell);
  }

  void tioga_set_highorder_callback_(void (*f1)(int*, int*),
                                     void (*f2)(int *,int *,double *),
                                     void (*f3)(int *,double *,int *,double *),
                                     void (*f4)(int *,double *,int *,int *,double *,double *,int *),
                                     void (*f5)(int *,int *,double *,int *,int *,double *))
  {
    tg->setcallback(f1,f2,f3,f4,f5);
    //getNodesPerCell=f1;
    //getReceptorNodes=f2;
    //donorInclusionTest=f3;
    //donorWeights=f4;
    //convert_to_modal=f5;
  }

  void tioga_set_solver_callback(solver* _solver)
  {
    tg->setcallback(_solver);
  }


  void tioga_delete_(void)
  {
    delete [] tg;
    if (nc) free(nc);
    if (nv) free(nv);
    if (vconn) free(vconn);
  }
}
