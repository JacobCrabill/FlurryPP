#include "MeshBlock.h"

#define ROW 0
#define COLUMN 1
#define NFRAC 1331
void MeshBlock::getCellIblanks(void)
{
  int i;
  int n,nvert,m;
  int icell;
  int inode[8];
  int ncount,flag;

  icell=0;
  for(n=0;n<ntypes;n++)
    {
      nvert=nv[n];
      for(i=0;i<nc[n];i++)
	{
	  flag=1;
	  iblank_cell[icell]=1;
	  ncount=0;
	  for(m=0;m<nvert && flag;m++)
	    {
	      inode[m]=vconn[n][nvert*i+m]-BASE;
	      if (iblank[inode[m]]==0) 
		{
		  iblank_cell[icell]=0;
		  flag=0;
		}
	      ncount=ncount+(iblank[inode[m]]==-1);
	    }
	  
	  if (flag) 
	    {
	      if (ncount ==nvert) iblank_cell[icell]=-1;
	    }
	  icell++;
	}
    }
}

void MeshBlock::getInternalNodes(void)
{
  int i,m;
  //
  nreceptorCells=0;
  //
  if (ctag!=NULL) free(ctag);
  ctag=(int *)malloc(sizeof(int)*ncells);
  //
  for(i=0;i<ncells;i++)
      if (iblank_cell[i]==-1) ctag[nreceptorCells++]=i+1;
  //
  if (pointsPerCell!=NULL) free(pointsPerCell);
  pointsPerCell=(int *)malloc(sizeof(int)*nreceptorCells);
  //
  maxPointsPerCell=0;
  //
  for(i=0;i<nreceptorCells;i++)
    {
      get_nodes_per_cell(&(ctag[i]),&(pointsPerCell[i]));
      ntotalPoints+=pointsPerCell[i];
      maxPointsPerCell=max(maxPointsPerCell,pointsPerCell[i]);
    }
  //
  if (rxyz !=NULL) free(rxyz);
  //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
  rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);
  //
  m=0;
  for(i=0;i<nreceptorCells;i++)
    {
      get_receptor_nodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
      m+=(3*pointsPerCell[i]);
    }
}

void MeshBlock::getExtraQueryPoints(OBB *obc,
			       int *nints,int **intData,
			       int *nreals, double **realData)
{
  int i,j,k;
  int i3;
  double xd[3];
  int *inode;
  int iptr;
  int m;

  inode=(int *)malloc(sizeof(int)*ntotalPoints);
  *nints=*nreals=0; 
  for(i=0;i<ntotalPoints;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(rxyz[i3+k]-obc->xc[k])*obc->vec[j][k];

      if (fabs(xd[0]) <= obc->dxc[0] &&
	  fabs(xd[1]) <= obc->dxc[1] &&
	  fabs(xd[2]) <= obc->dxc[2])
	{
	  inode[*nints]=i;
	  (*nints)++;
	  (*nreals)+=3;

	}
    }
  //
  (*intData)=(int *)malloc(sizeof(int)*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  //
  m=0;
  for(i=0;i<*nints;i++)
    {
      i3=3*inode[i];
      (*intData)[i]=inode[i];
      (*realData)[m++]=rxyz[i3];
      (*realData)[m++]=rxyz[i3+1];
      (*realData)[m++]=rxyz[i3+2];
    }
  //
  free(inode);
}  

void MeshBlock::processPointDonors(void)
{
  int i,j,m;
  double *frac;
  int icell;
  int ndim;
  int ierr;
  //
  ndim=NFRAC;
  frac=(double *) malloc(sizeof(double)*ndim);
  ninterp2=0;
  //
  for(i=0;i<nsearch;i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]]==1) ninterp2++;
  //
  if (interpList2) {
    for(i=0;i<ninterp2;i++)
      {
	free(interpList2[i].inode);
	free(interpList2[i].weights);
      }
    free(interpList2);
  } 
  interpList2=(INTERPLIST *)malloc(sizeof(INTERPLIST)*ninterp2);
  //  
  //printf("nsearch=%d %d\n",nsearch,myid);
  m=0;
  for(i=0;i<nsearch;i++)
    {
      if (donorId[i] > -1 && iblank_cell[donorId[i]]==1) 
	{
	  icell=donorId[i]+BASE;
	  interpList2[m].inode=(int *) malloc(sizeof(int));		  
          interpList2[m].nweights=0;
	  donor_frac(&(icell),
		     &(xsearch[3*i]),
		     &(interpList2[m].nweights),
		     &(interpList2[m].inode[0]),
		     frac,
		     &(rst[3*i]),&ndim);
	  interpList2[m].weights=(double *)malloc(sizeof(double)*interpList2[m].nweights);
	  for(j=0;j<interpList2[m].nweights;j++)
	    interpList2[m].weights[j]=frac[j];
	  interpList2[m].receptorInfo[0]=isearch[2*i];
	  interpList2[m].receptorInfo[1]=isearch[2*i+1];
	  m++;
	}
    }
  free(frac);
}

void MeshBlock::getInterpolatedSolutionAtPoints(int *nints,int *nreals,int **intData,
						double **realData,
						double *q,
						int nvar, int interptype)
{
  int i;
  int k,m,inode;
  double weight;
  double *qq;
  int icount,dcount;
  //
  qq=(double *)malloc(sizeof(double)*nvar);
  //
  (*nints)=ninterp2;
  (*nreals)=ninterp2*nvar;
  if ((*nints)==0) return;
  //
  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  icount=dcount=0;
  //
  if (interptype==ROW)
    {    
      for(i=0;i<ninterp2;i++)
	{
	  for(k=0;k<nvar;k++) qq[k]=0;
	  inode=interpList2[i].inode[0]-BASE;
	  for(m=0;m<interpList2[i].nweights;m++)
	    {
	      weight=interpList2[i].weights[m];
	      //if (weight < 0 || weight > 1.0) {
	      //	traced(weight);
	//	printf("warning: weights are not convex\n");
	  //    }
	      for(k=0;k<nvar;k++)
		qq[k]+=q[inode+m*nvar+k]*weight;
	    }
	  (*intData)[icount++]=interpList2[i].receptorInfo[0];
	  (*intData)[icount++]=interpList2[i].receptorInfo[1];
	  for(k=0;k<nvar;k++)
	    (*realData)[dcount++]=qq[k];
	}
    }
  //
  // no column-wise storage for high-order data
  //
}
	
void MeshBlock::updatePointData(double *q,double *qtmp,int nvar,int interptype)  
{
  int i,j,k,n,m;
  double *qout;
  int index_out;
  int npts;
  //
  npts=NFRAC;
  qout=(double *)malloc(sizeof(double)*nvar*npts);	
  //
  m=0;
  for(i=0;i<nreceptorCells;i++)
    {
      convert_to_modal(&(ctag[i]),&(pointsPerCell[i]),&(qtmp[m]),&npts,
		       &index_out,qout);
      index_out-=BASE;
      k=0;
      for(j=0;j<npts;j++)
	for(n=0;n<nvar;n++)
	  {
	    q[index_out+j*nvar+n]=qout[k];
	    k++;
	  }
      m+=(pointsPerCell[i]*nvar);
     }
  free(qout);
}
