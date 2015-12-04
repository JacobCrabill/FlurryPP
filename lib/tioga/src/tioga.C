#include "tioga.h"

#include "MeshBlock.h"
#include "utils.h"

tioga::tioga() {
  mb = NULL;
  holeMap=NULL;
  pc=NULL;
  sendCount=NULL;
  recvCount=NULL;
  obblist=NULL;
  isym=2;
  ihigh=0;
}

void tioga::setResolutions(double *nres,double *cres) {
  mb->setResolutions(nres,cres);
}

/**
 * set communicator
 * and initialize a few variables
 */
void tioga::setCommunicator(MPI_Comm communicator, int id_proc, int nprocs)
{
  scomm=communicator;
  myid=id_proc;
  numprocs=nprocs;
  sendCount=(int *) malloc(sizeof(int)*numprocs);
  recvCount=(int *) malloc(sizeof(int)*numprocs);

  // only one mesh block per process for now
  // this can be changed at a later date
  // but will be a fairly invasive change
  nblocks=1;
  mb=new MeshBlock[1];

  // instantiate the parallel communication class
  pc=new parallelComm[1];
  pc->myid=myid;
  pc->scomm=scomm;
  pc->numprocs=numprocs;
}

/**
 * register grid data for each mesh block
 */
void tioga::registerGridData(int btag,int nnodes,double *xyz,int *iblank, int nwbc,int nobc,
			     int *wbcnode,int *obcnode,int ntypes, int *nv, int *nc, int **vconn)
{
  mb->setData(btag,nnodes,xyz,iblank,nwbc,nobc,wbcnode,obcnode,ntypes,nv,nc,vconn);
  mb->myid=myid;
  mytag=btag;
}

void tioga::profile(void)
{
  mb->preprocess();
  //mb->writeGridFile(myid);
  //mb->writeOBB(myid);
  //if (myid==4) mb->writeOutput(myid);
  //if (myid==4) mb->writeOBB(myid);
}

void tioga::performConnectivity(void)
{
  getHoleMap();
  exchangeBoxes();
  exchangeSearchData();
  mb->search();
  exchangeDonors();
  if (ihigh) {
    mb->getCellIblanks();
    mb->writeCellFile(myid);
  }
  //mb->writeOutput(myid);
  //tracei(myid);
}

void tioga::performConnectivityHighOrder(void)
{
  //cout << "getInternalNodes" << endl;
  mb->getInternalNodes();
  MPI_Barrier(MPI_COMM_WORLD);
  //cout << "exchangePointSearchData" << endl;
  exchangePointSearchData();
  mb->ihigh=1;
  //cout << "search" << endl;
  mb->search();
  mb->processPointDonors();
  iorphanPrint=1;
}

int tioga::findPointDonor(double *x_pt)
{
  return mb->findPointDonor(x_pt);
}

std::unordered_set<int> tioga::findCellDonors(double *bbox)
{
  return mb->findCellDonors(bbox);
}

void tioga::dataUpdate(int nvar,double *q,int interptype)
{
  int i,j,k,m;
  int nints;
  int nreals;
  int *integerRecords;
  double *realRecords;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;
  int *icount,*dcount;

  // initialize send and recv packets
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);

  pc->initPackets(sndPack,rcvPack);

  // get the interpolated solution now
  integerRecords=NULL;
  realRecords=NULL;
  mb->getInterpolatedSolution(&nints,&nreals,&integerRecords,&realRecords,q,nvar,interptype);

  // populate the packets
  for(i=0;i<nints;i++)
    {
      k=integerRecords[2*i];
      sndPack[k].nints++;
      sndPack[k].nreals+=nvar;
    }

  for(k=0;k<nsend;k++)
    {
     sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
     sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
     icount[k]=dcount[k]=0;
    }

  m=0;
  for(i=0;i<nints;i++)
    {
      k=integerRecords[2*i];
      sndPack[k].intData[icount[k]++]=integerRecords[2*i+1];
      for(j=0;j<nvar;j++)
	sndPack[k].realData[dcount[k]++]=realRecords[m++];
    }

  // communicate the data across
  pc->sendRecvPackets(sndPack,rcvPack);

  // decode the packets and update the data
  for(k=0;k<nrecv;k++)
  {
    m=0;
    for(i=0;i<rcvPack[k].nints;i++)
    {
      mb->updateSolnData(rcvPack[k].intData[i],&rcvPack[k].realData[m],q,nvar,interptype);
      m+=nvar;
    }
  }

  // release all memory
  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);
}

void tioga::writeData(int nvar,double *q,int interptype)
{
  //mb->writeGridFile(myid);
  mb->writeFlowFile(myid,q,nvar,interptype);
}

void tioga::writeCellFile(void)
{
  mb->writeCellFile(myid);
}

void tioga::getDonorCount(int *dcount,int *fcount)
{
  mb->getDonorCount(dcount,fcount);
}

void tioga::getDonorInfo(int *receptors,int *indices,double *frac,int *dcount)
{
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  int i;
  mb->getDonorInfo(receptors,indices,frac);
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);

  // map receptor indices to actual processor id here
  for(i=0;i<3*(*dcount);i+=3)
    receptors[i]=sndMap[receptors[i]];
}

tioga::~tioga()
{
  int i;
  if (mb) delete[] mb;
  if (holeMap)
    {
      for(i=0;i<nmesh;i++)
	if (holeMap[i].existWall) free(holeMap[i].sam);
      delete [] holeMap;
    }
  if (pc) delete[] pc;
  if (sendCount) free(sendCount);
  if (recvCount) free(recvCount);
  if (obblist) free(obblist);
  if (myid==0) printf("#tioga :successfully cleared all the memory accessed\n");
}

void tioga::dataUpdate_highorder(int nvar,double *q,int interptype)
{
  int i,j,k,m;
  int nints;
  int nreals;
  int *integerRecords;
  double *realRecords;
  int nsend,nrecv;
  int *sndMap,*rcvMap;
  PACKET *sndPack,*rcvPack;
  int *icount,*dcount;
  double *qtmp;
  int norphanPoint;
  int *itmp;
  FILE *fp;
  char ofname[10];

  // initialize send and recv packets
  fp=NULL;
  icount=dcount=NULL;
  integerRecords=NULL;
  realRecords=NULL;

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend==0) return;
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nrecv);

  pc->initPackets(sndPack,rcvPack);

  // get the interpolated solution now
  integerRecords=NULL;
  realRecords=NULL;
  mb->getInterpolatedSolutionAtPoints(&nints,&nreals,&integerRecords,
				      &realRecords,q,nvar,interptype);

  // populate the packets
  for(i=0;i<nints;i++)
  {
    k=integerRecords[2*i];
    sndPack[k].nints++;
    sndPack[k].nreals+=nvar;
  }

  for(k=0;k<nsend;k++)
  {
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
    icount[k]=dcount[k]=0;
  }

  m=0;
  for(i=0;i<nints;i++)
  {
    k=integerRecords[2*i];
    sndPack[k].intData[icount[k]++]=integerRecords[2*i+1];
    for(j=0;j<nvar;j++)
      sndPack[k].realData[dcount[k]++]=realRecords[m++];
  }

  // communicate the data across
  pc->sendRecvPackets(sndPack,rcvPack);

  // decode the packets and update the data
  qtmp=(double *)malloc(sizeof(double)*nvar*mb->ntotalPoints);
  itmp=(int *) malloc(sizeof(int)*mb->ntotalPoints);
  for(i=0;i<mb->ntotalPoints;i++) itmp[i]=0;

  for(k=0;k<nrecv;k++)
  {
    m=0;
    for(i=0;i<rcvPack[k].nints;i++)
    {
      for(j=0;j<nvar;j++)
      {
        itmp[rcvPack[k].intData[i]]=1;
        qtmp[rcvPack[k].intData[i]*nvar+j]=rcvPack[k].realData[m];
        m++;
      }
    }
  }

  norphanPoint=0;
  for(i=0;i<mb->ntotalPoints;i++)
  {
    if (itmp[i]==0) {
      if (fp==NULL)
      {
        sprintf(ofname,"orphan%d.dat",myid);
        fp=fopen(ofname,"w");
      }
      mb->outputOrphan(fp,i);
      norphanPoint++;
    }
  }

  if (fp!=NULL) fclose(fp);

  if (norphanPoint > 0 && iorphanPrint) {
    printf("Warning::number of orphans in %d = %d of %d\n",myid,norphanPoint,
           mb->ntotalPoints);
    iorphanPrint=0;
  }

  mb->updatePointData(q,qtmp,nvar,interptype);

  // release all memory
  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
  free(qtmp);
  free(itmp);
  if (integerRecords) free(integerRecords);
  if (realRecords) free(realRecords);
  if (icount) free(icount);
  if (dcount) free(dcount);
}


void tioga::exchangeBoxes(void)
{
  int i,j,k,m;
  int *alltags;
  int *sndMap;
  int *rcvMap;
  int nsend;
  int nrecv;
  PACKET *sndPack,*rcvPack;
  //
  alltags=(int *)malloc(sizeof(int)*numprocs);
  MPI_Allgather(&mytag, 1, MPI_INT, alltags,1,MPI_INT,MPI_COMM_WORLD);
  //
  // count number of other processors to communicate to
  // in overset grid scenario, usually you do not communicate
  // to mesh blocks that carry the same mesh tag (i.e. you don't
  // talk to your sister partitions)
  //
  nsend=nrecv=0;
  for(i=0;i<numprocs;i++) if (alltags[i] != mytag) nsend++;
  //
  // In general we communicate forward
  // and backward, separate lists are maintained for
  // flexibility
  //
  nrecv=nsend;
  sndMap=(int *)malloc(sizeof(int)*nsend);
  rcvMap=(int *)malloc(sizeof(int)*nrecv);
  //
  for(i=0,m=0;i<numprocs;i++)
    if (alltags[i]!=mytag)
      {
        sndMap[m]=rcvMap[m]=i;
  m++;
      }
  //
  pc->setMap(nsend,nrecv,sndMap,rcvMap);
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  //
  for(k=0;k<nsend;k++)
    {
      sndPack[k].nints=0;
      sndPack[k].nreals=15;
      sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
      m=0;
      for(i=0;i<3;i++)
  for(j=0;j<3;j++)
    sndPack[k].realData[m++]=mb->obb->vec[i][j];
      for(i=0;i<3;i++)
  sndPack[k].realData[m++]=mb->obb->xc[i];
      for(i=0;i<3;i++)
  sndPack[k].realData[m++]=mb->obb->dxc[i];
    }
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  if (obblist) free(obblist);
  obblist=(OBB *) malloc(sizeof(OBB)*nrecv);
  //
  for(k=0;k<nrecv;k++)
    {
      m=0;
      for(i=0;i<3;i++)
  for(j=0;j<3;j++)
    obblist[k].vec[i][j]=rcvPack[k].realData[m++];
      for(i=0;i<3;i++)
  obblist[k].xc[i]=rcvPack[k].realData[m++];
      for(i=0;i<3;i++)
  obblist[k].dxc[i]=rcvPack[k].realData[m++];
    }
  //
  m=0;
  for(k=0;k<nrecv;k++)
    {
      if ( obbIntersectCheck(mb->obb->vec,mb->obb->xc,mb->obb->dxc,
           obblist[k].vec,obblist[k].xc,obblist[k].dxc) ||
     obbIntersectCheck(obblist[k].vec,obblist[k].xc,obblist[k].dxc,
           mb->obb->vec,mb->obb->xc,mb->obb->dxc))
  {
    rcvMap[m]=sndMap[k];
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        obblist[m].vec[i][j]=obblist[k].vec[i][j];
    for(i=0;i<3;i++)
      obblist[m].xc[i]=obblist[k].xc[i];
    for(i=0;i<3;i++)
      obblist[m].dxc[i]=obblist[k].dxc[i];
    m++;
  }
    }
  nsend=nrecv=m;
  for(i=0;i<nsend;i++) sndMap[i]=rcvMap[i];
  //printf("%d %d %d\n",myid,nsend,nrecv);
  //
  pc->setMap(nsend,nrecv,sndMap,rcvMap);
  //
  // free local memory
  //
  free(alltags);
  free(sndMap);
  free(rcvMap);
  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);

}


void tioga::exchangeDonors(void)
{
  int i,j,k;
  int nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int meshtag,pointid,remoteid;
  double donorRes;
  double *receptorResolution;
  int *donorRecords;
  int ninterp;
  int nrecords;
  int *sndMap;
  int *rcvMap;
  int *icount;
  int *dcount;

  donorRecords=NULL;
  icount=NULL;
  dcount=NULL;
  receptorResolution=NULL;


  // get the processor map for sending
  // and receiving

  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  if (nsend == 0) return;
  icount=(int *)malloc(sizeof(icount)*nsend);
  dcount=(int *)malloc(sizeof(icount)*nsend);

  // create packets to send and receive
  // and initialize them to zero

  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);

  pc->initPackets(sndPack,rcvPack);

  // get the data to send now

  mb->getDonorPacket(sndPack,nsend);

  // exchange the data (comm1)

  pc->sendRecvPackets(sndPack,rcvPack);

  // packet received has all the donor data
  // use this to populate link lists per point

  mb->initializeDonorList();

  for(k=0;k<nrecv;k++)
  {
    int nints=0;
    for(i=0;i<rcvPack[k].nints/3;i++)
    {
      meshtag=rcvPack[k].intData[nints++];
      pointid=rcvPack[k].intData[nints++];
      remoteid=rcvPack[k].intData[nints++];
      donorRes=rcvPack[k].realData[i];
      mb->insertAndSort(pointid,k,meshtag,remoteid,donorRes);
    }
  }

  // figure out the state of each point now (i.e. if its a hole or fringe or field)

  mb->processDonors(holeMap,nmesh,&donorRecords,&receptorResolution,&nrecords);

  // free Send and recv data

  pc->clearPackets(sndPack,rcvPack);

  // count number of records in each
  // packet

  for(i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].nints++;
    sndPack[k].nreals++;
  }

  // allocate the data containers

  for(k=0;k<nsend;k++)
  {
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
    sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
    icount[k]=dcount[k]=0;
  }

  for (i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
    sndPack[k].realData[dcount[k]++]=receptorResolution[i];
  }

  // now communicate the data (comm2)

  pc->sendRecvPackets(sndPack,rcvPack);

  // create the interpolation list now

  ninterp=0;
  for(k=0;k<nrecv;k++)
    ninterp+=rcvPack[k].nints;

  mb->initializeInterpList(ninterp);

  ninterp=0;
  for(k=0;k<nrecv;k++)
    for(i=0;i<rcvPack[k].nints;i++)
      mb->findInterpData(&ninterp,rcvPack[k].intData[i],rcvPack[k].realData[i]);
  mb->set_ninterp(ninterp);

  //printf("process %d has (%d,%d) points to interpolate out %d donors\n",myid,ninterp,m,mb->donorCount);

  pc->clearPackets(sndPack,rcvPack);
  if (donorRecords) free(donorRecords);
  donorRecords=NULL;

  // cancel donors that have conflict

  mb->getCancellationData(&nrecords,&donorRecords);
  //printf("process %d has %d cancelled receptors\n",myid,nrecords);


  for(i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].nints++;
  }
  for(k=0;k<nsend;k++)
  {
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
    icount[k]=0;
  }
  for(i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
  }

  // communicate the cancellation data (comm 3)

  pc->sendRecvPackets(sndPack,rcvPack);

  for(k=0;k<nrecv;k++)
  {
    for(i=0;i<rcvPack[k].nints;i++)
      mb->cancelDonor(rcvPack[k].intData[i]);
  }

  if (donorRecords) free(donorRecords);
  donorRecords=NULL;
  pc->clearPackets(sndPack,rcvPack);

  mb->getInterpData(&nrecords,&donorRecords);

  for(i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].nints++;
  }
  for(k=0;k<nsend;k++)
  {
    sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
    icount[k]=0;
  }
  for(i=0;i<nrecords;i++)
  {
    k=donorRecords[2*i];
    sndPack[k].intData[icount[k]++]=donorRecords[2*i+1];
  }

  // communicate the final receptor data (comm 4)

  pc->sendRecvPackets(sndPack,rcvPack);
  mb->clearIblanks();
  for(k=0;k<nrecv;k++)
    for(i=0;i<rcvPack[k].nints;i++)
      mb->setIblanks(rcvPack[k].intData[i]);

  // finished the communication, free all
  // memory now

  if (donorRecords) free(donorRecords);
  if (receptorResolution) free(receptorResolution);
  if (icount) free(icount);
  if (dcount) free(dcount);

  // free Send and recv data

  pc->clearPackets(sndPack,rcvPack);
  free(sndPack);
  free(rcvPack);
}


void tioga::exchangeSearchData(void)
{
  int i,j,k,l,m;
  int icount,dcount;
  int nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int *sndMap;
  int *rcvMap;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);

  for(i=0;i<nsend;i++)
    {
      sndPack[i].nints=sndPack[i].nreals=0;
      sndPack[i].intData=NULL;
      sndPack[i].realData=NULL;
    }

  for(i=0;i<nrecv;i++)
    {
      rcvPack[i].nints=rcvPack[i].nreals=0;
      rcvPack[i].intData=NULL;
      rcvPack[i].realData=NULL;
    }
  //
  // now get data for each packet
  //
  for(k=0;k<nsend;k++)
    mb->getQueryPoints(&obblist[k],
           &sndPack[k].nints,&sndPack[k].intData,
           &sndPack[k].nreals,&sndPack[k].realData);
  //
  // exchange the data
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // now assort the data into the search
  // list arrays
  //
  mb->nsearch=0;
  for(k=0;k<nrecv;k++)
    mb->nsearch+=rcvPack[k].nints;
  //
  // if these were already allocated
  // get rid of them
  //
  if (mb->xsearch) free(mb->xsearch);
  if (mb->isearch) free(mb->isearch);
  if (mb->donorId) free(mb->donorId);
  //
  // allocate query point storage
  //
  mb->xsearch=(double *)malloc(sizeof(double)*3*mb->nsearch);
  mb->isearch=(int *)malloc(2*sizeof(int)*mb->nsearch);
  mb->donorId=(int *)malloc(sizeof(int)*mb->nsearch);
  //
  // now fill the query point arrays
  //
  icount=0;
  dcount=0;
  for(k=0;k<nrecv;k++)
  {
    l=0;
    for(j=0;j<rcvPack[k].nints;j++)
    {
      mb->isearch[icount++]=k;
      mb->isearch[icount++]=rcvPack[k].intData[j];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
    }
  }
  for(i=0;i<nsend;i++)
  {
    if (sndPack[i].nints > 0) free(sndPack[i].intData);
    if (sndPack[i].nreals >0) free(sndPack[i].realData);
  }
  for(i=0;i<nrecv;i++)
  {
    if (rcvPack[i].nints > 0) free(rcvPack[i].intData);
    if (rcvPack[i].nreals >0) free(rcvPack[i].realData);
  }
  free(sndPack);
  free(rcvPack);
  //printf("%d %d\n",myid,mb->nsearch);
}

//
// routine for extra query points
// have to unify both routines here
// FIX later ...
//
void tioga::exchangePointSearchData(void)
{
  int i,j,k,l,m;
  int icount,dcount;
  int nsend,nrecv;
  PACKET *sndPack,*rcvPack;
  int *sndMap;
  int *rcvMap;
  //
  // get the processor map for sending
  // and receiving
  //
  pc->getMap(&nsend,&nrecv,&sndMap,&rcvMap);
  //
  // create packets to send and receive
  // and initialize them to zero
  //
  sndPack=(PACKET *)malloc(sizeof(PACKET)*nsend);
  rcvPack=(PACKET *)malloc(sizeof(PACKET)*nrecv);
  //
  for(i=0;i<nsend;i++)
    {
      sndPack[i].nints=sndPack[i].nreals=0;
      sndPack[i].intData=NULL;
      sndPack[i].realData=NULL;
    }
  //
  for(i=0;i<nrecv;i++)
    {
      rcvPack[i].nints=rcvPack[i].nreals=0;
      rcvPack[i].intData=NULL;
      rcvPack[i].realData=NULL;
    }
  //
  // now get data for each packet
  //
  for(k=0;k<nsend;k++)
    mb->getExtraQueryPoints(&obblist[k],
          &sndPack[k].nints,&sndPack[k].intData,
          &sndPack[k].nreals,&sndPack[k].realData);
  MPI_Barrier(MPI_COMM_WORLD);
  //if (myid==0) printf("AAAAA\n");
  //
  // exchange the data
  //
  pc->sendRecvPackets(sndPack,rcvPack);
  //
  // now assort the data into the search
  // list arrays
  //
  mb->nsearch=0;
  for(k=0;k<nrecv;k++)
    mb->nsearch+=rcvPack[k].nints;
  //
  // if these were already allocated
  // get rid of them
  //
  if (mb->xsearch) free(mb->xsearch);
  if (mb->isearch) free(mb->isearch);
  if (mb->donorId) free(mb->donorId);
  if (mb->rst) free(mb->rst);
  //
  // allocate query point storage
  //
  mb->xsearch=(double *)malloc(sizeof(double)*3*mb->nsearch);
  mb->isearch=(int *)malloc(2*sizeof(int)*mb->nsearch);
  mb->donorId=(int *)malloc(sizeof(int)*mb->nsearch);
  mb->rst=(double *) malloc(sizeof(double)*3*mb->nsearch);
  //
  // now fill the query point arrays
  //
  icount=0;
  dcount=0;
  for(k=0;k<nrecv;k++)
  {
    l=0;
    for(j=0;j<rcvPack[k].nints;j++)
    {
      mb->isearch[icount++]=k;
      mb->isearch[icount++]=rcvPack[k].intData[j];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
      mb->xsearch[dcount++]=rcvPack[k].realData[l++];
    }
  }
  for(i=0;i<nsend;i++)
  {
    if (sndPack[i].nints > 0) free(sndPack[i].intData);
    if (sndPack[i].nreals >0) free(sndPack[i].realData);
  }
  for(i=0;i<nrecv;i++)
  {
    if (rcvPack[i].nints > 0) free(rcvPack[i].intData);
    if (rcvPack[i].nreals >0) free(rcvPack[i].realData);
  }
  free(sndPack);
  free(rcvPack);
  //printf("%d %d\n",myid,mb->nsearch);
}


/**
 * Create hole maps for all grids
 * this routine is not efficient
 * since it does mutiple global reduce ops
 * have to change it at a later date when
 * there is more time to develop code
 */
void tioga::getHoleMap(void)
{
  int i,j,k,m;
  int ii,jj,kk;
  double wbox[6];
  int existWall;
  int meshtag,maxtag;
  int *existHoleLocal;
  int *existHole;
  double *bboxLocal;
  double *bboxGlobal;
  double ds[3],dsmax,dsbox;
  int bufferSize;
  FILE *fp;
  char fname[80];
  char intstring[7];
 //
 // get the local bounding box
 //
 mb->getWallBounds(&meshtag,&existWall,wbox);
 MPI_Allreduce(&meshtag,&maxtag,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 //
 if (holeMap)
 {
   for(i=0;i<nmesh;i++)
     if (holeMap[i].existWall) free(holeMap[i].sam);
   delete [] holeMap;
 }
 holeMap=new HOLEMAP[maxtag];
 //
 existHoleLocal=(int *)malloc(sizeof(int)*maxtag);
 existHole=(int *)malloc(sizeof(int)*maxtag);
 //
 for(i=0;i<maxtag;i++) existHole[i]=existHoleLocal[i]=0;
 //
 existHoleLocal[meshtag-1]=existWall;
 //
 MPI_Allreduce(existHoleLocal,existHole,maxtag,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 //
 for(i=0;i<maxtag;i++) holeMap[i].existWall=existHole[i];
 //
 bboxLocal=(double *) malloc(sizeof(double)*6*maxtag);
 bboxGlobal=(double *) malloc(sizeof(double)*6*maxtag);
 //
 for(i=0;i<3*maxtag;i++) bboxLocal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxLocal[i+3*maxtag]=-BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i+3*maxtag]=-BIGVALUE;

 //
 for(i=0;i<3;i++)
   {
     bboxLocal[3*(meshtag-1)+i]=wbox[i];
     bboxLocal[3*(meshtag-1)+i+3*maxtag]=wbox[i+3];
   }
 //
 // get the global bounding box info across all the
 // partitions for all meshes
 //
 MPI_Allreduce(bboxLocal,bboxGlobal,3*maxtag,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 MPI_Allreduce(&(bboxLocal[3*maxtag]),&(bboxGlobal[3*maxtag]),3*maxtag,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 //
 // find the bounding box for each mesh
 // from the globally reduced data
 //
 for(i=0;i<maxtag;i++)
 {
   if (holeMap[i].existWall)
   {
     for(j=0;j<3;j++)
     {
       holeMap[i].extents[j]=bboxGlobal[3*i+j];
       holeMap[i].extents[j+3]=bboxGlobal[3*i+j+3*maxtag];
       ds[j]=holeMap[i].extents[j+3]-holeMap[i].extents[j];
     }
     dsmax=max(ds[0],ds[1]);
     dsmax=max(dsmax,ds[2]);
     dsbox=dsmax/64;

     for(j=0;j<3;j++)
     {
       holeMap[i].extents[j]-=(2*dsbox);
       holeMap[i].extents[j+3]+=(2*dsbox);
       holeMap[i].nx[j]=floor((double)max((holeMap[i].extents[j+3]-holeMap[i].extents[j])/dsbox,1.));
     }
     bufferSize=holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
     holeMap[i].sam=(int *)malloc(sizeof(int)*bufferSize);
     holeMap[i].samLocal=(int *)malloc(sizeof(int)*bufferSize);
     for(j=0;j<bufferSize;j++) holeMap[i].sam[j]=holeMap[i].samLocal[j]=0;
   }
 }
 //
 // mark the wall boundary cells in the holeMap
 //
 if (holeMap[meshtag-1].existWall) {
   mb->markWallBoundary(holeMap[meshtag-1].samLocal,holeMap[meshtag-1].nx,holeMap[meshtag-1].extents);
 }
 //
 // allreduce the holeMap of each mesh
 //
 for(i=0;i<maxtag;i++)
 {
   if (holeMap[i].existWall)
   {
     bufferSize=holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
     MPI_Allreduce(holeMap[i].samLocal,holeMap[i].sam,bufferSize,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   }
 }
 //
 // set the global number of meshes to maxtag
 //
 nmesh=maxtag;
 //
 // now fill the holeMap
 //
 for(i=0;i<maxtag;i++)
   if (holeMap[i].existWall) fillHoleMap(holeMap[i].sam,holeMap[i].nx,isym);
 //
 // output the hole map
 //
 //this->outputHoleMap();
 //
 // free local memory
 //
 free(existHoleLocal);
 free(existHole);
 free(bboxLocal);
 free(bboxGlobal);
}

/**
 * Output the hole map to a tecplot compatible file
*/
void tioga::outputHoleMap(void)
{
  int i,k;
  int nnodes,ncells;
  int ns1,ns2;
  int ii,jj,kk,m;
  FILE *fp;
  double ds[3];
  char intstring[7];
  char fname[80];

  for(i=0;i<nmesh;i++)
    if (holeMap[i].existWall)
       {
   sprintf(intstring,"%d",100000+i+100*myid);
   sprintf(fname,"holeMap%s.dat",&(intstring[1]));
   fp=fopen(fname,"w");
   fprintf(fp,"TITLE =\"Tioga output\"\n");
   fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
   nnodes=(holeMap[i].nx[0]+1)*(holeMap[i].nx[1]+1)*(holeMap[i].nx[2]+1);
   ncells=(holeMap[i].nx[0])*(holeMap[i].nx[1])*(holeMap[i].nx[2]);
   fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,ncells);
   fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");
   for(k=0;k<3;k++) ds[k]=(holeMap[i].extents[k+3]-holeMap[i].extents[k])/(holeMap[i].nx[k]);
   //
   for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
     for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
       for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
         fprintf(fp,"%.14e\n",ii*ds[0]);
   for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
     for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
       for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
         fprintf(fp,"%.14e\n",jj*ds[1]);
   for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
     for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
       for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
         fprintf(fp,"%.14e\n",kk*ds[2]);
   m=0;
   for(kk=0;kk<holeMap[i].nx[2];kk++)
     for(jj=0;jj<holeMap[i].nx[1];jj++)
       for(ii=0;ii<holeMap[i].nx[0];ii++)
         {
     fprintf(fp,"%f\n",(double)holeMap[i].sam[m]);
     m++;
         }

   m=0;
         ns1=holeMap[i].nx[0]+1;
   ns2=(holeMap[i].nx[1]+1)*ns1;
   for(kk=0;kk<holeMap[i].nx[2];kk++)
     for(jj=0;jj<holeMap[i].nx[1];jj++)
       for(ii=0;ii<holeMap[i].nx[0];ii++)
         {
     m=kk*ns2+jj*ns1+ii+1;
     fprintf(fp,"%d %d %d %d %d %d %d %d\n",m,m+1,m+1+ns1,m+ns1,
       m+ns2,m+1+ns2,m+ns2+ns1+1,m+ns1+ns2);
         }
       }
 fclose(fp);
}

void tioga::set_cell_iblank(int *iblank_cell)
{
  mb->set_cell_iblank(iblank_cell);
}

//! Set callback functions for high-order processing. See MeshBlock.h for details
void tioga::setcallback(void (*f1)(int*, int*),
      void (*f2)(int *,int *,double *),
      void (*f3)(int *,double *,int *,double *),
      void (*f4)(int *,double *,int *,int *,double *,double *,int *),
      void (*f5)(int *,int *,double *,int *,int*,double *))
{
  mb->setcallback(f1,f2,f3,f4,f5);
  ihigh=1;
}

void tioga::setcallback(solver* _solver)
{
  mb->setcallback(_solver);
  ihigh=1;
}
