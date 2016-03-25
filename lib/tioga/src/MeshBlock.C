#include "MeshBlock.h"

#include "utils.h"
#include "funcs.hpp"

void MeshBlock::setData(int btag,int nnodesi,double *xyzi, int *ibli,int nwbci, int nobci,
			int *wbcnodei,int *obcnodei,
			int ntypesi,int *nvi,int *nci,int **vconni)
{
  // set internal pointers
  meshtag=btag;
  nnodes=nnodesi;
  x=xyzi;
  iblank=ibli;
  nwbc=nwbci;
  nobc=nobci;
  wbcnode=wbcnodei;
  obcnode=obcnodei;

  ntypes=ntypesi;

  nv=nvi;
  nc=nci;
  vconn=vconni;

  ncells=0;
  for (int i = 0; i < ntypes; i++) ncells += nc[i];
}

void MeshBlock::preprocess(void)
{
  // set all iblanks = 1
  for (int i = 0; i < nnodes; i++) iblank[i] = 1;

  // find oriented bounding boxes
  if (obb) free(obb);
  obb = (OBB *) malloc(sizeof(OBB));
  findOBB(x,obb->xc,obb->dxc,obb->vec,nnodes);
  tagBoundary();
}

void MeshBlock::tagBoundary(void)
{
  vector<int> inode;
  double xv[8][3];
  std::vector<int> iflag(nnodes, 0);

  // do this only once
  // i.e. when the meshblock is first
  // initialized, cellRes would be NULL in this case

  if (cellRes == NULL)
  {
    cellRes = (double *) malloc(sizeof(double)*ncells);
    nodeRes = (double *) malloc(sizeof(double)*nnodes);

    if (userSpecifiedNodeRes == NULL && userSpecifiedCellRes == NULL)
    {
      for (int i = 0; i < nnodes; i++) nodeRes[i] = 0.0;

      int k = 0;
      for (int n = 0; n < ntypes; n++)
      {
        int nvert = nv[n];
        inode.resize(nvert);
        for (int i = 0; i < nc[n]; i++)
        {
          for (int m = 0; m < nvert; m++)
          {
            inode[m] = vconn[n][nvert*i+m]-BASE;
            if (m < 8) {
              int i3 = 3*inode[m];
              for (int j = 0; j < 3; j++)
                xv[m][j] = x[i3+j];
            }
          }
          double vol = computeCellVolume(xv,min(nvert,8));  /// TODO: add support for high-order hexas
          cellRes[k++] = vol;
          for (int m = 0; m < nvert; m++)
          {
            iflag[inode[m]]++;
            nodeRes[inode[m]] += vol;
          }
        }
      }
    }
    else
    {
      int k = 0;
      for (int n = 0; n < ntypes; n++)
      {
        for (int i = 0; i < nc[n]; i++)
        {
          cellRes[k]=userSpecifiedCellRes[k];
          k++;
        }
      }
      for (int k = 0; k < nnodes; k++) nodeRes[k] = userSpecifiedNodeRes[k];
    }

    // compute nodal resolution as the average of
    // all the cells associated with it. This takes care
    // of partition boundaries as well.

    for (int i = 0; i < nnodes; i++)
    {
      if (iflag[i] != 0)  nodeRes[i] /= iflag[i];
      iflag[i] = 0;
    }

    // now tag the overset-boundary nodes
    // reuse the iflag array

    for (int i = 0; i < nobc; i++)
    {
      iflag[(obcnode[i]-BASE)] = 1;
    }

    // now tag all the nodes of overset-boundary cells
    // to be mandatory receptors

    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for (int i = 0; i < nc[n]; i++)
      {
        int itag = 0;
        for (int m = 0; m < nvert; m++)
        {
          inode[m] = vconn[n][nvert*i+m]-BASE;
          if (iflag[inode[m]]) itag = 1;
        }
        if (itag)
        {
          for (int m = 0; m < nvert; m++)
          {
            nodeRes[inode[m]]=BIGVALUE;
          }
        }
      }
    }

    // now tag all the cells which have
    // mandatory receptors as nodes as not acceptable
    // donors
    int k = 0;
    for (int n = 0; n < ntypes; n++)
    {
      int nvert = nv[n];
      for (int i = 0; i < nc[n]; i++)
      {
        for (int m = 0; m < nvert; m++)
        {
          inode[m] = vconn[n][nvert*i+m]-BASE;
          if (iflag[inode[m]])
          {
            cellRes[k] = BIGVALUE;
            break;
          }
        }
        k++;
      }
    }
  }
}

void MeshBlock::writeGridFile(int bid)
{
  char fname[80];
  char intstring[7];
  FILE *fp;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"part%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
	  ncells);
  for (int i = 0; i < nnodes; i++)
  {
    fprintf(fp,"%.14e %.14e %.14e %d\n",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
  }

  int ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    int nvert_type = min(nvert,8);
    for (int i = 0; i < nc[n]; i++)
    {
      if (nvert_type==4)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba);
      }
      else if (nvert_type==5)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba);
      }
      else if (nvert_type==6)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+5]+ba);
      }
      else if (nvert_type==8)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+6]+ba,
            vconn[n][nvert*i+7]+ba);
      }
    }
  }
  fclose(fp);
  return;
}

void MeshBlock::writeCellFile(int bid)
{
  char fname[80];
  char intstring[7];
  FILE *fp;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"cell%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\",\"IBLANK_CELL\" ");
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,
	  ncells);
  fprintf(fp,"VARLOCATION =  (1=NODAL, 2=NODAL, 3=NODAL, 4=NODAL,5=CELLCENTERED)\n");
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i+1]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%lf\n",x[3*i+2]);
  for (int i = 0; i < nnodes; i++) fprintf(fp,"%d\n",iblank[i]);
  for (int i = 0; i < ncells; i++) fprintf(fp,"%d\n",iblank_cell[i]);

  int ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    int nvert_type = min(nvert,8);
    for (int i = 0; i < nc[n]; i++)
    {
      if (nvert_type==4)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba);
      }
      else if (nvert_type==5)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba);
      }
      else if (nvert_type==6)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+5]+ba);
      }
      else if (nvert_type==8)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+6]+ba,
            vconn[n][nvert*i+7]+ba);
      }
    }
  }
  fclose(fp);
  return;
}

void MeshBlock::writeFlowFile(int bid,double *q,int nvar,int type)
{
  char fname[80];
  char qstr[2];
  char intstring[7];
  FILE *fp;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"flow%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Tioga output\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\" ");
  for (int i = 0; i < nvar; i++)
  {
    sprintf(qstr,"Q%d",i);
    fprintf(fp,"\"%s\",",qstr);
  }
  fprintf(fp,"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",nnodes,
          ncells);

  if (type==0)
  {
    for (int i = 0; i < nnodes; i++)
    {
      fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
      for (int j = 0; j < nvar; j++)
        fprintf(fp,"%lf ",q[i*nvar+j]);
      //for (int j = 0; j < nvar; j++)
      //  fprintf(fp,"%lf ", x[3*i]+x[3*i+1]+x[3*i+2]);
      fprintf(fp,"\n");
    }
  }
  else
  {
    for (int i = 0; i < nnodes; i++)
    {
      fprintf(fp,"%lf %lf %lf %d ",x[3*i],x[3*i+1],x[3*i+2],iblank[i]);
      for (int j = 0; j < nvar; j++)
        fprintf(fp,"%lf ",q[j*nnodes+i]);
      fprintf(fp,"\n");
    }
  }

  int ba=1-BASE;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert = nv[n];
    int nvert_type = min(nvert,8);
    for (int i = 0; i < nc[n]; i++)
    {
      if (nvert_type==4)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+3]+ba);
      }
      else if (nvert_type==5)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+4]+ba);
      }
      else if (nvert_type==6)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+5]+ba);
      }
      else if (nvert_type==8)
      {
        fprintf(fp,"%d %d %d %d %d %d %d %d\n",
                vconn[n][nvert*i]+ba,
            vconn[n][nvert*i+1]+ba,
            vconn[n][nvert*i+2]+ba,
            vconn[n][nvert*i+3]+ba,
            vconn[n][nvert*i+4]+ba,
            vconn[n][nvert*i+5]+ba,
            vconn[n][nvert*i+6]+ba,
            vconn[n][nvert*i+7]+ba);
      }
    }
  }
  fclose(fp);
  return;
}

void MeshBlock::getWallBounds(int &mtag,int *existWall, double wbox[6])
{
  int i,j,i3;
  int inode;

  mtag = meshtag+(1-BASE);
  if (nwbc <=0) {
    *existWall=0;
    for (int i = 0; i < 6; i++) wbox[i]=0;
    return;
  }

  *existWall=1;
  wbox[0]=wbox[1]=wbox[2]=BIGVALUE;
  wbox[3]=wbox[4]=wbox[5]=-BIGVALUE;

  for (int i = 0; i < nwbc; i++)
    {
      inode=wbcnode[i]-BASE;
      i3=3*inode;
      for (int j = 0; j < 3; j++)
	{
	  wbox[j]=min(wbox[j],x[i3+j]);
	  wbox[j+3]=max(wbox[j+3],x[i3+j]);
	}
    }

}

void MeshBlock::markWallBoundary(int *sam,int nx[3],double extents[6])
{
  vector<int> iflag(ncells,0);
  vector<int> inode(nnodes,0);
  double ds[3];
  int imin[3];
  int imax[3];

  for (int i = 0; i<nwbc; i++)
  {
    int ii=wbcnode[i]-BASE;
    inode[ii]=1;
  }

  // mark wall boundary cells
  int m = 0;
  for (int n = 0; n<ntypes; n++)
  {
    int nvert=nv[n];
    for (int i = 0; i<nc[n]; i++)
    {
      for (int j = 0; j<nvert; j++)
      {
        int ii=vconn[n][nvert*i+j]-BASE;
        if (inode[ii]==1)
        {
          iflag[m]=1;
          break;
        }
      }
      m++;
    }
  }

  // find delta's in each directions
  for (int k = 0; k<3; k++) ds[k]=(extents[k+3]-extents[k])/nx[k];

  // mark sam cells with wall boundary cells now
  m = 0;
  for (int n = 0; n<ntypes; n++)
  {
    int nvert = nv[n];
    for (int i = 0; i<nc[n]; i++)
    {
      if (iflag[m]==1)
      {
        // find the index bounds of each wall boundary cell
        // bounding box
        imin[0]=imin[1]=imin[2]=BIGINT;
        imax[0]=imax[1]=imax[2]=-BIGINT;
        for (int j = 0; j < nvert; j++)
        {
          int i3=3*(vconn[n][nvert*i+j]-BASE);
          for (int k = 0; k < 3; k++)
          {
            double xv=x[i3+k];
            int iv=floor((xv-extents[k])/ds[k]);
            imin[k]=min(imin[k],iv);
            imax[k]=max(imax[k],iv);
          }
        }
        for (int j = 0; j < 3; j++)
        {
          imin[j]=max(imin[j],0);
          imax[j]=min(imax[j],nx[j]-1);
        }

        // mark sam to 1
        for (int kk = imin[2]; kk < imax[2]+1; kk++)
          for (int jj = imin[1]; jj < imax[1]+1; jj++)
            for (int ii = imin[0]; ii < imax[0]+1; ii++)
            {
              int mm = kk*nx[1]*nx[0]+jj*nx[0]+ii;
              sam[mm] = 2;
            }
      }
      m++;
    }
  }
}

void MeshBlock::getQueryPoints(OBB *obc,int &nints,int **intData,int &nreals, double **realData)
{
  double xd[3];
  vector<int> inode(nnodes);

  nints = nreals = 0;
  for (int i = 0; i < nnodes; i++)
  {
    int i3=3*i;
    for (int j = 0; j < 3; j++) xd[j]=0;
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        xd[j]+=(x[i3+k]-obc->xc[k])*obc->vec[j][k];

    if (fabs(xd[0]) <= obc->dxc[0] && fabs(xd[1]) <= obc->dxc[1] && fabs(xd[2]) <= obc->dxc[2])
    {
      inode[nints]=i;
      nints++;
      nreals += 3;

    }
  }

  (*intData)=(int *)malloc(sizeof(int)*nints);
  (*realData)=(double *)malloc(sizeof(double)*nreals);

  int m = 0;
  for (int i = 0; i < nints; i++)
  {
    int i3=3*inode[i];
    (*intData)[i]=inode[i];
    (*realData)[m++]=x[i3];
    (*realData)[m++]=x[i3+1];
    (*realData)[m++]=x[i3+2];
  }
}

void MeshBlock::writeOBB(int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"box%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,1);

  for (int l = 0; l < 2; l++)
  {
    il=2*(l%2)-1;
    for (int k = 0; k < 2; k++)
    {
      ik=2*(k%2)-1;
      for (int j = 0; j < 2; j++)
      {
        ij=2*(j%2)-1;
        xx[0]=xx[1]=xx[2]=0;
        for (int m = 0; m < 3; m++)
          xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
              +ik*obb->vec[1][m]*obb->dxc[1]
              +il*obb->vec[2][m]*obb->dxc[2];
        fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
      }
    }
  }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fprintf(fp,"%e %e %e\n",obb->xc[0],obb->xc[1],obb->xc[2]);
  for (int k = 0; k < 3; k++)
   fprintf(fp,"%e %e %e\n",obb->vec[0][k],obb->vec[1][k],obb->vec[2][k]);
  fprintf(fp,"%e %e %e\n",obb->dxc[0],obb->dxc[1],obb->dxc[2]);
  fclose(fp);
}

/*! destructor that deallocates all the
 *  the dynamic objects inside
 */
MeshBlock::~MeshBlock()
{
  int i;
  //
  // free all data that is owned by this MeshBlock
  // i.e not the pointers of the external code.
  //
  if (cellRes) free(cellRes);
  if (nodeRes) free(nodeRes);
  if (elementBbox) free(elementBbox);
  if (elementList) free(elementList);
  if (adt) delete[] adt;
  if (donorList) {
    for (int i = 0; i < nnodes; i++) deallocateLinkList(donorList[i]);
    free(donorList);
  }
  if (interpList) {
    for (int i = 0; i < ninterp; i++)
      {
	free(interpList[i].inode);
	free(interpList[i].weights);
      }
    free(interpList);
  }
  if (interpList2) {
    for (int i = 0; i < ninterp2; i++)
      {
	free(interpList2[i].inode);
	free(interpList2[i].weights);
      }
    free(interpList2);
  }
  if (obb) free(obb);
  if (isearch) free(isearch);
  if (xsearch) free(xsearch);
  if (rst) free(rst);
  if (interp2donor) free(interp2donor);
  if (cancelList) deallocateLinkList2(cancelList);
  if (ctag) free(ctag);
  if (pointsPerCell) free(pointsPerCell);
  if (rxyz) free(rxyz);
  // need to add code here for other objects as and
  // when they become part of MeshBlock object
};
//
// set user specified node and cell resolutions
//
void MeshBlock::setResolutions(double *nres,double *cres)
{
  userSpecifiedNodeRes=nres;
  userSpecifiedCellRes=cres;
}

int MeshBlock::getCellIndex(int adtEle)
{
  return elementList[adtEle];
}


void MeshBlock::checkContainment(int *cellIndex, int adtElement, double *xsearch)
{
  double xv[8][3];
  double frac[8];

  int icell = elementList[adtElement];

  // get the type of the cell & the type-specific cell index
  int isum = 0;
  int ic = -1;
  int N = 0;
  for (int n = 0; n < ntypes; n++)
  {
    isum += nc[n];
    if (icell < isum)
    {
      ic = icell-(isum-nc[n]);
      N = n;
      break;
    }
  }

  if (ic<0) {
    std::cout << "invalid icell (adtElement) in checkContainment" << std::endl;
    exit(1);
  }

  *cellIndex = icell+BASE;

  if (ihigh==0)
  {
    /* --- Normal [Not high-order] --- */

    int nvert = nv[N];
    if (nvert > 8)
    {
      /* --- Quadratic or higher-order shape functions - use general func --- */
      //Solver->donorInclusionTest(cellIndex,xsearch,&isInEle,&refloc[0]);

      vector<double> xv2(nvert*3);
      for (int m = 0; m < nvert; m++)
      {
        int i3 = 3*(vconn[N][nvert*ic+m]-BASE);
        for (int j = 0; j < 3; j++)
          xv2[m*3+j] = x[i3+j];
      }

      double refloc[3];
      bool isInEle = getRefLocNewton(xv2.data(), xsearch, &refloc[0], nvert, 3);

      // if any of the nodal weights are not in between [-TOL 1+TOL] discard cell
      if (!isInEle)
      {
        *cellIndex=-1;
        return;
      }
      return;
    }
    else
    {
      /* --- Linear shape functions --- */
      // now collect all the vertices in the array xv
      for (int m = 0; m < nvert; m++)
      {
        int i3 = 3*(vconn[N][nvert*ic+m]-BASE);
        for (int j = 0; j < 3; j++)
          xv[m][j] = x[i3+j];
      }

      computeNodalWeights(xv,xsearch,frac,nvert);

      // if any of the nodal weights are not in between [-TOL 1+TOL] discard cell
      for (int m = 0; m < nvert; m++)
      {
        if ((frac[m]+TOL)*(frac[m]-1.0-TOL) > 0)
        {
          *cellIndex=-1;
          return;
        }
      }
      return;
    }
  }
  else
  {
    /* --- Quadratic or higher-order shape functions - use general func --- */

    int nvert = nv[N];
    vector<double> xv2(nvert*3);
    for (int m = 0; m < nvert; m++)
    {
      int i3 = 3*(vconn[N][nvert*ic+m]-BASE);
      for (int j = 0; j < 3; j++)
        xv2[m*3+j] = x[i3+j];
    }

    double refloc[3];
    bool isInEle = getRefLocNewton(xv2.data(), xsearch, &refloc[0], nvert, 3);

    if (!isInEle) *cellIndex = -1;

//    /* --- High-Order: Use user-given callback function --- */
//    int icell1=icell+BASE;
//    *cellIndex=-1;
//    bool isInEle = false;
//    Solver->donorInclusionTest(&icell1,xsearch,&isInEle,&(rst[ipoint]));
//    if (!passFlag) *cellIndex = -1;
    return;
  }

}

void MeshBlock::getInterpolatedSolution(int* nints, int* nreals, int** intData, double** realData,
                                        double* q, int nvar, int interptype)
{
  int i;
  int k,m,inode;
  double weight;
  double* qq;
  int icount,dcount;

  qq=(double *)malloc(sizeof(double)*nvar);

  (*nints)=(*nreals)=0;
  for (int i = 0; i<ninterp; i++)
  {
    if (!interpList[i].cancel)
    {
      (*nints)++;
      (*nreals)=(*nreals)+nvar;
    }
  }
  if ((*nints)==0) return;

  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  icount=dcount=0;

  if (interptype==ROW) // Row-major vs. col-major data
  {
    for (int i = 0; i<ninterp; i++)
    {
      if (!interpList[i].cancel)
      {
        for (int k = 0; k<nvar; k++) qq[k]=0;
        for (int m = 0; m<interpList[i].nweights; m++)
        {
          inode=interpList[i].inode[m];
          weight=interpList[i].weights[m];
          if (weight < 0 || weight > 1.0) {
            traced(weight);
            printf("warning: weights are not convex\n");
          }
          for (int k = 0; k<nvar; k++)
            qq[k]+=q[inode*nvar+k]*weight;
        }
        (*intData)[icount++]=interpList[i].receptorInfo[0];
        (*intData)[icount++]=interpList[i].receptorInfo[1];
        for (int k = 0; k<nvar; k++)
          (*realData)[dcount++]=qq[k];
      }
    }
  }
  else if (interptype==COLUMN)
  {
    for (int i = 0; i<ninterp; i++)
    {
      if (!interpList[i].cancel)
      {
        for (int k = 0; k<nvar; k++) qq[k]=0;
        for (int m = 0; m<interpList[i].nweights; m++)
        {
          inode=interpList[i].inode[m];
          weight=interpList[i].weights[m];
          for (int k = 0; k<nvar; k++)
            qq[k]+=q[k*nnodes+inode]*weight;
        }
        (*intData)[icount++]=interpList[i].receptorInfo[0];
        (*intData)[icount++]=interpList[i].receptorInfo[1];
        for (int k = 0; k<nvar; k++)
          (*realData)[dcount++]=qq[k];
      }
    }
  }
}

void MeshBlock::updateSolnData(int inode,double* qvar,double* q,int nvar,int interptype)
{
  int k;

  if (interptype==ROW) // Row-major vs. col-major data
  {
    for (int k = 0; k<nvar; k++)
      q[inode* nvar+k]=qvar[k];
  }
  if (interptype==COLUMN)
  {
    for (int k = 0; k<nvar; k++)
      q[nnodes* k+inode]=qvar[k];
  }
}

void MeshBlock::getDonorCount(int* dcount,int* fcount)
{
  int i;
  *dcount=0;
  *fcount=0;
  for (int i = 0; i<ninterp; i++)
  {
    if (!interpList[i].cancel)
    {
      (*dcount)++;
      (*fcount)+=interpList[i].nweights;
    }
  }
}

void MeshBlock::getDonorInfo(int* receptors,int* indices,double* frac)
{
  int i,j,k,m;
  int dcount=0;

  j = 0;
  k = 0;
  for (int i = 0; i<ninterp; i++)
  {
    if (!interpList[i].cancel)
    {
      for (int m = 0; m<interpList[i].nweights; m++)
      {
        indices[j]=interpList[i].inode[m];
        frac[j]=interpList[i].weights[m];
        j++;
      }
      receptors[k++]=interpList[i].receptorInfo[0];
      receptors[k++]=interpList[i].receptorInfo[1];
      receptors[k++]=interpList[i].nweights;
    }
  }
}


/*! Use nodal iblank values to determine iblank value for whole cell */
void MeshBlock::getCellIblanks(void)
{
  int i;
  int n,nvert,m;
  int icell;
  int inode[8];
  int ncount,flag;

  icell=0;
  for (int n = 0; n < ntypes; n++)
  {
    nvert=nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      flag=1;
      iblank_cell[icell] = NORMAL;
      ncount=0;
      for (int m = 0; m < nvert && flag; m++)
      {
        // If any node is set to hole (iblank == 0), set entire cell to hole
        inode[m]=vconn[n][nvert*i+m]-BASE;
        if (iblank[inode[m]] == HOLE)
        {
          iblank_cell[icell] = HOLE;
          flag=0;
        }
        if (iblank[inode[m]] == FRINGE) ncount++;
      }

      // If all nodes of cell are fringe nodes, set cell to fringe (-1)
      if (flag)
      {
        if (ncount == nvert) iblank_cell[icell]= FRINGE;
      }
      icell++;
    }
  }
}

void MeshBlock::getInternalNodes(void)
{
  nreceptorCells=0;

  if (ctag!=NULL) free(ctag);
  ctag=(int *)malloc(sizeof(int)*ncells);

  for (int i = 0; i < ncells; i++)
    if (iblank_cell[i]==-1) ctag[nreceptorCells++]=i+BASE;

  if (pointsPerCell!=NULL) free(pointsPerCell);
  pointsPerCell=(int *)malloc(sizeof(int)*nreceptorCells);

  maxPointsPerCell=0;

  for (int i = 0; i < nreceptorCells; i++)
  {
    // Use user-specified callback function to assign number of points within cell
    Solver->getNodesPerCell(&(ctag[i]),&(pointsPerCell[i]));
    ntotalPoints+=pointsPerCell[i];
    maxPointsPerCell=max(maxPointsPerCell,pointsPerCell[i]);
  }

  if (rxyz !=NULL) free(rxyz);
  //printf("getInternalNodes : %d %d\n",myid,ntotalPoints);
  rxyz=(double *)malloc(sizeof(double)*ntotalPoints*3);

  int m = 0;
  for (int i = 0; i < nreceptorCells; i++)
  {
    Solver->getReceptorNodes(&(ctag[i]),&(pointsPerCell[i]),&(rxyz[m]));
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
  int m;

  inode=(int *)malloc(sizeof(int)*ntotalPoints);
  *nints=*nreals=0;
  for (int i = 0; i < ntotalPoints; i++)
  {
    i3=3*i;
    for (int j = 0; j < 3; j++) xd[j]=0;

    // Calculate distance of receptor point from center of OBB, rotated to align
    // with the OBB's axes
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        xd[j]+=(rxyz[i3+k]-obc->xc[k])*obc->vec[j][k];

    // If receptor point is within bounding box, add to list
    if (fabs(xd[0]) <= obc->dxc[0] && fabs(xd[1]) <= obc->dxc[1] && fabs(xd[2]) <= obc->dxc[2])
    {
      inode[*nints]=i;
      (*nints)++;
      (*nreals)+=3;

    }
  }

  (*intData)=(int *)malloc(sizeof(int)*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));

  m = 0;
  for (int i = 0; i < *nints; i++)
  {
    i3=3*inode[i];
    (*intData)[i]=inode[i];
    (*realData)[m++]=rxyz[i3];
    (*realData)[m++]=rxyz[i3+1];
    (*realData)[m++]=rxyz[i3+2];
  }

  free(inode);
}

void MeshBlock::processPointDonors(void)
{
  double *frac;
  int    *ind;
  int icell;
  int ndim;

  ndim=NFRAC;
  frac=(double *) malloc(sizeof(double)*ndim);
  ind =(int    *) malloc(sizeof(int   )*ndim);
  ninterp2=0;

  for (int i = 0; i<nsearch; i++)
    if (donorId[i] > -1 && iblank_cell[donorId[i]]==NORMAL) ninterp2++;

  if (interpList2) {
    for (int i = 0; i<ninterp2; i++)
    {
      free(interpList2[i].inode);
      free(interpList2[i].weights);
    }
    free(interpList2);
  }
  interpList2=(INTERPLIST *)malloc(sizeof(INTERPLIST)*ninterp2);

  int m = 0;
  for (int i = 0; i<nsearch; i++)
  {
    if (donorId[i] > -1 && iblank_cell[donorId[i]]==1)
    {
      icell=donorId[i]+BASE;
      interpList2[m].inode=(int *) malloc(sizeof(int));
      interpList2[m].nweights=0;

      // Use the user-specified callback function to get the donor cell's interpolation weights
      // Outputs: nweights, inode, frac,
      Solver->donorWeights(&(icell), &(xsearch[3*i]), &(interpList2[m].nweights), ind, frac, &(rst[3*i]), &ndim);

      interpList2[m].weights=(double *)malloc(sizeof(double)*interpList2[m].nweights);
      interpList2[m].inode  =(int    *)malloc(sizeof(int   )*interpList2[m].nweights);

      for (int j = 0; j<interpList2[m].nweights; j++) {
        interpList2[m].weights[j] = frac[j];
        interpList2[m].inode[j]   = ind[j];
      }
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

  qq=(double *)malloc(sizeof(double)*nvar);

  (*nints)=ninterp2;
  (*nreals)=ninterp2*nvar;
  if ((*nints)==0) return;

  (*intData)=(int *)malloc(sizeof(int)*2*(*nints));
  (*realData)=(double *)malloc(sizeof(double)*(*nreals));
  icount=dcount=0;

  if (interptype==ROW)
  {
    for (int i = 0; i < ninterp2; i++)
    {
      for (int k = 0; k < nvar; k++) qq[k]=0;
      inode=interpList2[i].inode[0]-BASE;
      for (int m = 0; m < interpList2[i].nweights; m++)
      {
        weight=interpList2[i].weights[m];
        //if (weight < 0 || weight > 1.0) {
        //	traced(weight);
        //	printf("warning: weights are not convex\n");
        //    }
        for (int k = 0; k < nvar; k++)
          qq[k]+=q[inode+m*nvar+k]*weight;
      }
      (*intData)[icount++]=interpList2[i].receptorInfo[0];
      (*intData)[icount++]=interpList2[i].receptorInfo[1];
      for (int k = 0; k < nvar; k++)
        (*realData)[dcount++]=qq[k];
    }
  }
  //
  // no column-wise storage for high-order data
  //
}

void MeshBlock::updatePointData(double *q,double *qtmp,int nvar,int interptype)
{
  int index_out;

  int npts = NFRAC;
  double *qout = (double *)malloc(sizeof(double)*nvar*npts);

  int m = 0;
  for (int i = 0; i < nreceptorCells; i++)
  {
    Solver->convertToModal(&(ctag[i]),&(pointsPerCell[i]),&(qtmp[m]),&npts,&index_out,qout);
    index_out-=BASE;
    int k = 0;
    for (int j = 0; j < npts; j++) {
      for (int n = 0; n < nvar; n++)
      {
        q[index_out+j*nvar+n]=qout[k];
        k++;
      }
    }
    m+=(pointsPerCell[i]*nvar);
  }

  free(qout);
}

/* ---- Donor-Search Functions ---- */

void MeshBlock::search(void)
{
  double xd[3] = {0,0,0};
  double dxc[3] = {0,0,0};
  double xmin[3] = {0,0,0};
  double xmax[3] = {0,0,0};

  // form the bounding box of the query points
  /* --- For Flurry: Must have access to the ADT for all ranks, so don't return! --- */
//  if (nsearch == 0) {
//    donorCount=0;
//    return;
//  }

  OBB* obq=(OBB *) malloc(sizeof(OBB));
  findOBB(xsearch,obq->xc,obq->dxc,obq->vec,nsearch);

  // find all the cells that may have intersections with
  // the OBB

  int* icell=(int *)malloc(sizeof(int)*ncells);
  for (int i = 0; i < ncells; i++) icell[i]=-1;

  int iptr=-1;
  int cell_count=0;
  for (int n = 0; n < ntypes; n++)
  {
    int nvert=nv[n];
    for (int i = 0; i < nc[n]; i++)
    {
      // find each cell that has
      // overlap with the bounding box
      xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
      xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
      for (int m = 0; m < nvert; m++)
      {
        int i3=3*(vconn[n][nvert*i+m]-BASE);
        for (int j = 0; j < 3; j++)
        {
          xd[j]=0;
          for (int k = 0; k < 3; k++)
            xd[j]+=(x[i3+k]-obq->xc[k])*obq->vec[j][k];
          xmin[j]=min(xmin[j],xd[j]);
          xmax[j]=max(xmax[j],xd[j]);
        }
        for (int j = 0; j < 3; j++)
        {
          xd[j]=(xmax[j]+xmin[j])*0.5;
          dxc[j]=(xmax[j]-xmin[j])*0.5;
        }
      }
      if (fabs(xd[0]) <= (dxc[0]+obq->dxc[0]) &&
          fabs(xd[1]) <= (dxc[1]+obq->dxc[1]) &&
          fabs(xd[2]) <= (dxc[2]+obq->dxc[2]))
      {

        /* Create a LIFO stack ('icell') with all the cells
         * that have a bounding-box intersection with the
         * query-point bounding box
         */
        icell[i]=iptr;
        iptr=i;
        cell_count++;
      }
    }
  }

  // now find the axis aligned bounding box
  // of each cell in the LIFO stack to build the
  // ADT

  if (elementBbox) free(elementBbox);
  if (elementList) free(elementList);
  elementBbox=(double *)malloc(sizeof(double)*cell_count*6);
  elementList=(int *)malloc(sizeof(int)*cell_count);

  int k=iptr;
  int l=0;
  int p=0;
  while(k!=-1)
  {
    int cellindex=k;
    int isum = 0;
    int ic = 0;
    int N = 0;
    for (int n = 0; n < ntypes; n++)
    {
      isum+=nc[n];
      if (cellindex < isum)
      {
        ic = cellindex-(isum-nc[n]);
        N = n;
        break;
      }
    }
    int nvert=nv[N];
    xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
    xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
    for (int m = 0; m < nvert; m++)
    {
      int i3=3*(vconn[N][nvert*ic+m]-BASE);
      for (int j = 0; j < 3; j++)
      {
        xmin[j]=min(xmin[j],x[i3+j]);
        xmax[j]=max(xmax[j],x[i3+j]);
      }
    }

    elementBbox[l++]=xmin[0];
    elementBbox[l++]=xmin[1];
    elementBbox[l++]=xmin[2];
    elementBbox[l++]=xmax[0];
    elementBbox[l++]=xmax[1];
    elementBbox[l++]=xmax[2];

    elementList[p++]=k;

    k=icell[k];
  }

  // build the ADT now

  if (adt) {
    adt->clearData();
  } else {
    adt=new ADT[1];
  }

  int ndim=6;
  adt->buildADT(ndim,cell_count,elementBbox);

  if (donorId) free(donorId);
  donorId=(int*)malloc(sizeof(int)*nsearch);

  donorCount=0;
  ipoint=0;
  for (int i = 0; i < nsearch; i++)
  {
    adt->searchADT_point(this,&(donorId[i]),&(xsearch[3*i]));
    if (donorId[i] > -1) {
      donorCount++;
    }
    ipoint+=3;
  }

  free(icell);
  free(obq);
}

int MeshBlock::findPointDonor(double *x_pt)
{
  int foundCell;
  adt->searchADT_point(this,&foundCell,x_pt);
  return foundCell;
}

std::unordered_set<int> MeshBlock::findCellDonors(double *bbox)
{
  std::unordered_set<int> foundCells;
  adt->searchADT_box(elementList,foundCells,bbox);
  return foundCells;
}

/* ---- Bookkeeping Functions ---- */


void MeshBlock::getDonorPacket(PACKET *sndPack, int nsend)
{
  int i,k;
  int *icount;
  int *dcount;

  icount=(int *)malloc(sizeof(int)*nsend);
  dcount=(int *)malloc(sizeof(int)*nsend);
  //
  // count numbers to send first
  //
  for (int i = 0; i < nsearch; i++)
    {
      if (donorId[i] > -1)
        {
          k=isearch[2*i];
          sndPack[k].nints+=3;
          sndPack[k].nreals++;
        }
    }
  for (int k = 0; k < nsend; k++)
    {
      sndPack[k].intData=(int *)malloc(sizeof(int)*sndPack[k].nints);
      sndPack[k].realData=(double *)malloc(sizeof(double)*sndPack[k].nreals);
    }

  for (int i = 0; i < nsend; i++) {icount[i]=dcount[i]=0;};
  //
  for (int i = 0; i < nsearch; i++)
    {
      if (donorId[i] > -1)
        {
          k=isearch[2*i];
          sndPack[k].intData[icount[k]++]=meshtag;               // mesh tag
          sndPack[k].intData[icount[k]++]=isearch[2*i+1];        // point id
          sndPack[k].intData[icount[k]++]=i;                         // point id on the donor side
          sndPack[k].realData[dcount[k]++]=cellRes[donorId[i]];  // donor resolution
        }
    }
  free(icount);
  free(dcount);
}
void MeshBlock::initializeDonorList(void)
{
  int i;
  if (donorList)
    {
      for (int i = 0; i < nnodes; i++)
  deallocateLinkList(donorList[i]);
      free(donorList);
    }
  donorList=(DONORLIST **)malloc(sizeof(DONORLIST *)*nnodes);
  for (int i = 0; i < nnodes; i++)
    donorList[i]=NULL;
}

void MeshBlock::insertAndSort(int pointid,int senderid,int meshtagdonor, int remoteid,
            double donorRes)
{
  DONORLIST *temp1;
  temp1=(DONORLIST *)malloc(sizeof(DONORLIST));
  temp1->donorData[0]=senderid;
  temp1->donorData[1]=meshtagdonor;
  temp1->donorData[2]=remoteid;
  temp1->donorRes=donorRes;
  insertInList(&donorList[pointid],temp1);
}

void MeshBlock::processDonors(HOLEMAP *holemap, int nmesh, int **donorRecords,double **receptorResolution,
            int *nrecords)
{
  DONORLIST *temp;
  int meshtagdonor;

  // first mark hole points

  int* iflag=(int *)malloc(sizeof(int)*nmesh);

  for (int i = 0; i < nnodes; i++)
  {
    iblank[i] = NORMAL;
    if (donorList[i]==NULL)
    {
      for (int j = 0; j < nmesh; j++) {
        if (j!=(meshtag-BASE) && holemap[j].existWall)
        {
          if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam.data(),holemap[j].extents))
          {
            iblank[i] = HOLE;
            break;
          }
        }
      }
    }
    else
    {
      temp=donorList[i];
      for (int j = 0; j < nmesh; j++) iflag[j]=0;
      while(temp!=NULL)
      {
        meshtagdonor=temp->donorData[1]-BASE;
        iflag[meshtagdonor]=1;
        temp=temp->next;
      }
      for (int j = 0; j < nmesh; j++)
      {
        if (j!=(meshtag-BASE) && holemap[j].existWall)
        {
          if (!iflag[j]) {
            if (checkHoleMap(&x[3*i],holemap[j].nx,holemap[j].sam.data(),holemap[j].extents))
            {
              iblank[i] = HOLE;
              break;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < nwbc; i++) {
    if (iblank[wbcnode[i]-BASE]==HOLE) {
      printf("--------------------------------------------------------------------\n");
      printf("Alarm from process %d : wall node is being tagged as a hole %d %p\n",myid,wbcnode[i]-BASE,
             donorList[wbcnode[i]-BASE]);
      int ii=wbcnode[i]-BASE;
      printf("xloc=%e %e %e\n",x[3*ii],x[3*ii+1],x[3*ii+2]);
      printf("Computations will continue, but may suffer from accuracy problems\n");
      printf("Please recheck positions of your grids\n");
      printf("--------------------------------------------------------------------\n");
    }
  }

  // mark mandatory fringes as neighbors (up to nfringe depth)
  // of hole points

  int* mtag = new int[nnodes];
  int* mtag1 = new int[nnodes];

  for (int i = 0; i < nnodes; i++)
  {
    mtag[i]=mtag1[i]=0;
    if (iblank[i] == HOLE) mtag[i]=mtag1[i]=NORMAL;
  }

  for (int iter=0; iter<nfringe; iter++)
  {
    for (int n = 0; n < ntypes; n++)
    {
      int nvert=nv[n];
      for (int i = 0; i < nc[n]; i++)
      {
        for (int m = 0; m < nvert; m++)
        {
          if (mtag[(vconn[n][nvert*i+m]-BASE)]==NORMAL)
          {
            for (int mm = 0;mm<nvert;mm++)
              if (m!=mm && mtag[vconn[n][nvert*i+mm]-BASE] !=NORMAL)
                mtag1[vconn[n][nvert*i+mm]-BASE] = NORMAL;
          }
        }
      }
    }
    for (int i = 0; i < nnodes; i++) mtag[i]=mtag1[i];
  }

  for (int i = 0; i < nnodes; i++)
    if (mtag1[i] && iblank[i]) nodeRes[i]=BIGVALUE;

  delete [] mtag;
  delete [] mtag1;

  // now find fringes

  *nrecords=0;
  for (int i = 0; i < nnodes; i++)
  {
    if (donorList[i]!=NULL && iblank[i]!=HOLE)
    {
      temp=donorList[i];
      while(temp!=NULL)
      {
        if (temp->donorRes < nodeRes[i])
        {
          iblank[i] = FRINGE;
          (*nrecords)++;
          break;
        }
        temp=temp->next;
      }
    }
  }

  // set the records to send back to the donor
  // process

  (*donorRecords)=(int *)malloc(sizeof(int)*2*(*nrecords));
  (*receptorResolution)=(double *)malloc(sizeof(double)*(*nrecords));

  int m = 0;
  int k = 0;
  for (int i = 0; i < nnodes; i++)
  {
    if (iblank[i] == FRINGE)
    {
      temp=donorList[i];
      (*receptorResolution)[k++]=nodeRes[i];
      (*donorRecords)[m++]=temp->donorData[0];
      (*donorRecords)[m++]=temp->donorData[2];
    }
  }

  // release local memory
  free(iflag);
}


/* ---- Data Update Functions ---- */

void MeshBlock::initializeInterpList(int ninterp_input)
{
  int i;

  if (interpList) {
    for (int i = 0; i < ninterp; i++)
      {
  free(interpList[i].inode);
  free(interpList[i].weights);
      }
    free(interpList);
  }
  ninterp=ninterp_input;
  interpList=(INTERPLIST *)malloc(sizeof(INTERPLIST)*ninterp);
  if (cancelList) deallocateLinkList2(cancelList);
  cancelList=NULL;
  ncancel=0;
  if (interp2donor) free(interp2donor);
  interp2donor=(int *)malloc(sizeof(int)*nsearch);
  for (int i = 0; i < nsearch; i++) interp2donor[i]=-1;

}

void MeshBlock::findInterpData(int *recid,int irecord,double receptorRes)
{
  int idonor = 0;
  double xv[8][3];
  vector<double> xv2;
  double xp[3];
  double frac[8];
  vector<double> frac2;
  int inode[8];
  INTEGERLIST *clist;

  int procid=isearch[2*irecord];
  int pointid=isearch[2*irecord+1];
  int i3 = 3*irecord;
  xp[0]=xsearch[i3];
  xp[1]=xsearch[i3+1];
  xp[2]=xsearch[i3+2];

  // Get the donor cell for the given receptor point
  int isum = 0;
  int N = 0;
  for (int n = 0; n < ntypes; n++)
  {
    isum+=nc[n];
    if (donorId[irecord] < isum)
    {
      idonor=donorId[irecord]-(isum-nc[n]);
      N = n;
      break;
    }
  }

  // Get the physical node positions of the donor cell
  int acceptFlag = 1;
  int nvert = nv[N];

  if (nvert > 8)
    xv2.resize(nvert*3);

  for (int iv=0; iv<nvert; iv++)
  {
    inode[iv]=vconn[N][nvert*idonor+iv]-BASE;
    int i3=3*inode[iv];
    if (iblank[inode[iv]] <=0) {
      // If donor-cell node is also a receptor or hole node, skip
      acceptFlag=0;
    }

    if (nvert > 8)
      for (int j = 0; j < 3; j++)
        xv2[iv*3+j]=x[i3+j];
    else
      for (int j = 0; j < 3; j++)
        xv[iv][j]=x[i3+j];
  }

  if (acceptFlag==0 && receptorRes!=BIGVALUE) return;
  if (receptorRes==BIGVALUE)
  {
    clist=cancelList;

    // go to the end of the list
    if (clist != NULL) while(clist->next != NULL) clist=clist->next;

    for (int iv=0; iv<nvert; iv++)
    {
      inode[iv]=vconn[N][nvert*idonor+iv]-BASE;
      if (iblank[inode[iv]]==-1 && nodeRes[inode[iv]]!=BIGVALUE)
      {
        iblank[inode[iv]]=1;
        if (clist == NULL)
        {
          clist=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
          clist->inode=inode[iv];
          clist->next=NULL;
          cancelList=clist;
        }
        else
        {
          clist->next=(INTEGERLIST *)malloc(sizeof(INTEGERLIST));
          clist->next->inode=inode[iv];
          clist->next->next=NULL;
          clist=clist->next;
        }
        ncancel++;
      }
    }
  }

  /* Compute the (linear) shape-function interpolation weights of the receptor
   * node within the donor cell */
  if (nvert <= 8) {
    computeNodalWeights(xv,xp,frac,nvert); /// TODO: support nvert > 8
  } else {
    double xref[3];
    int isInEle;
    Solver->donorInclusionTest(&idonor,&xp[0],&isInEle,&xref[0]);

    frac2.resize(nvert);
    shape_hex(point(&xref[0]),frac2,nvert);
  }

  interp2donor[irecord]=*recid;
  interpList[*recid].cancel=0;
  interpList[*recid].nweights=nvert;
  interpList[*recid].receptorInfo[0]=procid;
  interpList[*recid].receptorInfo[1]=pointid;
  interpList[*recid].inode=(int *)malloc(sizeof(int)*nvert);
  interpList[*recid].weights=(double *)malloc(sizeof(double)*nvert);
  for (int iv=0;iv<nvert;iv++)
  {
    interpList[*recid].inode[iv]=inode[iv];
    if (nvert > 8)
      interpList[*recid].weights[iv]=frac2[iv];
    else
      interpList[*recid].weights[iv]=frac[iv];
  }
  (*recid)++;
}

void MeshBlock::set_ninterp(int ninterp_input)
{
  ninterp=ninterp_input;
}

void MeshBlock::getCancellationData(int *nrecords,int **intData)
{
  int i;
  int inode;
  INTEGERLIST *clist;
  *nrecords=ncancel;
  if (ncancel > 0)
    {
      (*intData)=(int *)malloc(sizeof(int)*(*nrecords)*2);
      i = 0;
      for (clist=cancelList;clist!=NULL;clist=clist->next)
  {
    inode=clist->inode;
    (*intData)[i++]=donorList[inode]->donorData[0];
    (*intData)[i++]=donorList[inode]->donorData[2];
  }
    }
}

void MeshBlock::cancelDonor(int irecord)
{
  int iptr;
  iptr=interp2donor[irecord];
  if (iptr > -1) interpList[iptr].cancel=1;
}

void MeshBlock::getInterpData(int *nrecords, int **intData)
{
  int i,k;
  //
  *nrecords=0;
  for (int i = 0; i < ninterp; i++)
    if (!interpList[i].cancel) (*nrecords)++;
  //
  (*intData)=(int *)malloc(sizeof(int)*2*(*nrecords));
  for (i=0,k = 0; i < ninterp; i++)
    if (!interpList[i].cancel) {
       (*intData)[k++]=interpList[i].receptorInfo[0];
       (*intData)[k++]=interpList[i].receptorInfo[1];
    }
}

void MeshBlock::clearIblanks(void)
{
  int i;
  for (int i = 0; i < nnodes; i++)
     if (iblank[i] < 0) iblank[i]=1;
}

void MeshBlock::setIblanks(int inode)
{
/*  if (fabs(nodeRes[inode]-BIGVALUE) < TOL)
    {
      iblank[inode]=-2;
    }
  else*/
//    {
   iblank[inode]=-1;
//    }
}
