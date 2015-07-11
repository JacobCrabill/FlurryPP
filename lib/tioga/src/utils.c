#include "utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "codetypes.h"

extern "C" {
void cellvolume_(double* vol, double xv[8][3], int numverts[4][6], int faceInfo[4][24], int* nfaces, int* nvert);
}

extern void kaiser_wrap_(double *,int *,int *,double *,double *,double *,int *);
/***
 ** find oriented bounding box for a given set of points
*/
void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes)
{
  int i,j,k,m,i3;
  double *aa;
  double *eigenv;
  double trace,sume;
  int ier;
  double xd[3];
  double xmin[3],xmax[3];
  int nrows,ncols;
  //
  xc[0]=xc[1]=xc[2]=0;
  //
  // find centroid coordinates (not true centroid)
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      xc[0]+=x[i3];
      xc[1]+=x[i3+1];
      xc[2]+=x[i3+2];
    }
  //
  xc[0]/=nnodes;
  xc[1]/=nnodes;
  xc[2]/=nnodes;
  if (nnodes <4) {
    vec[0][0]=vec[1][1]=vec[2][2]=1;
    vec[0][1]=vec[1][0]=vec[1][2]=vec[2][1]=0;
    vec[0][2]=vec[2][0]=0;

    if (nnodes==1) {
	    dxc[0]=1e-3;
      dxc[1]=1e-3;
      dxc[2]=1e-3;
      return;
    }
    else if (nnodes==2)	{
      dxc[0]=max(1e-3,fabs(x[3]-x[0]))*0.5;
      dxc[1]=max(1e-3,fabs(x[4]-x[1]))*0.5;
      dxc[2]=max(1e-3,fabs(x[5]-x[2]))*0.5;
      return;
    }
    else
    {
      for(i=0;i<nnodes;i++) {
        i3=3*i;
        for(j=0;j<3;j++)
          dxc[j]=max(1e-3,fabs(x[i3+j]-x[0]));
      }
      return;
    }
  }
  //
  // find co-variance matrix
  // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
  //
  aa=(double *) malloc(sizeof(double)*9);
  eigenv=(double *)malloc(sizeof(double)*3);
  //
  for(i=0;i<9;i++) aa[i]=0;
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      aa[0]+=((x[i3]-xc[0])*(x[i3]-xc[0]));
      aa[4]+=((x[i3+1]-xc[1])*(x[i3+1]-xc[1]));
      aa[8]+=((x[i3+2]-xc[2])*(x[i3+2]-xc[2]));
      aa[3]+=((x[i3]-xc[0])*(x[i3+1]-xc[1]));
      aa[6]+=((x[i3]-xc[0])*(x[i3+2]-xc[2]));
      aa[7]+=((x[i3+1]-xc[1])*(x[i3+2]-xc[2]));
    }
  aa[1]=aa[3];
  aa[2]=aa[6];
  aa[5]=aa[7];
  //
  // use kaisers method to estimate
  // eigen values and vectors of the covariance matrix
  //
  nrows=3;
  ncols=3;
  kaiser_wrap_(aa,&nrows,&ncols,eigenv,&trace,&sume,&ier);
  //
  // copy the eigen vector basis on to vec
  //
  m=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
       vec[i][j]=aa[m++];
      }
  //
  // find min and max bounds in the bounding box
  // vector basis
  //
  for(j=0;j<3;j++)
    {
      xmax[j]=-BIGVALUE;
      xmin[j]=BIGVALUE;
    }
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      //
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(x[i3+k]-xc[k])*vec[j][k];
      //
      for(j=0;j<3;j++)
	{
	  xmax[j]=max(xmax[j],xd[j]);
	  xmin[j]=min(xmin[j],xd[j]);
	}
    }
  //
  // find the extents of the box
  // and coordinates of the center w.r.t. xc
  // increase extents by 1% for tolerance
  //
  for(j=0;j<3;j++)
    {
      dxc[j]=(xmax[j]-xmin[j])*0.5*1.01;
      xd[j]=(xmax[j]+xmin[j])*0.5;
    }
  //
  // find the center of the box in
  // actual cartesian coordinates
  //
  for(j=0;j<3;j++)
    {
    for(k=0;k<3;k++)
      xc[j]+=(xd[k]*vec[k][j]);
    }
  //
  free(aa);
  free(eigenv);
}
/**
 check if a point is inside the
 provided hole map
*/
int checkHoleMap(double *x,int *nx,int *sam,double *extents)
{
  int i;
  int mm;
  double dx[3];
  int ix[3];

  for(i=0;i<3;i++) dx[i]=(extents[i+3]-extents[i])/nx[i];
  for(i=0;i<3;i++)
    {
      ix[i]=(x[i]-extents[i])/dx[i];
      if (ix[i] < 0 || ix[i] > nx[i]-1) return 0;
    }
  mm=ix[2]*nx[1]*nx[0]+ix[1]*nx[0]+ix[0];
  return sam[mm];
}

/**
 fill a given hole map using iterative
 flood fill from outside the marked boundary.
 boundary is marked by "2"
*/
void fillHoleMap(int *holeMap, int ix[3],int isym)
{
  int m;
  int ii,jj,kk,mm;
  int ns2;
  int i,j,k;
  int ipaint,npaint,nneig;
  int mk[6];
  //
  // now start from outside and paint the
  // the exterior
  //
  ns2=ix[0]*ix[1];
  //
  for(kk=0;kk<ix[2];kk+=(ix[2]-1))
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii++)
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj+=(ix[1]-1))
      for(ii=0;ii<ix[0];ii++)
        {
  	  mm=kk*ns2+jj*ix[0]+ii;
         holeMap[mm]=1;
  	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii+=(ix[0]-1))
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}
  npaint=ns2*ix[2];
  while(npaint > 0)
    {
      npaint=0;
      for(k=1;k<ix[2]-1;k++)
        for(j=1;j<ix[1]-1;j++)
          for(i=1;i<ix[0]-1;i++)
            {
              m=k*ns2+j*ix[0]+i;
              if (holeMap[m]==0)
                {
                  ipaint=0;
                  if (isym==1)
                   {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[2]=m-ix[0];
		     mk[3]=m+ix[0];
		     mk[4]=m-1;
		     mk[5]=m+1;
		     nneig=4;
		   }
		  else if (isym==2)
		    {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[4]=m-ix[0];
		     mk[5]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else if (isym==3)
		    {
		     mk[4]=m-ns2;
		     mk[5]=m+ns2;
		     mk[0]=m-ix[0];
		     mk[1]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else
		    {
		      mk[0]=m-ns2;
		      mk[1]=m+ns2;
		      mk[2]=m-ix[0];
		      mk[3]=m+ix[0];
		      mk[4]=m-1;
		      mk[5]=m+1;
		      nneig=6;
		    }
                  for (kk=0;kk<nneig && ipaint==0;kk++)
		    {
		      ipaint=(ipaint || holeMap[mk[kk]]==1);
		    }
                  if (ipaint > 0)
                    {
                      holeMap[m]=1;
                      npaint++;
                    }
                }
            }
    }
  for(i=0;i<ix[2]*ix[1]*ix[0];i++)
   {
    holeMap[i]=(holeMap[i] ==0 || holeMap[i]==2);
   }
}

int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
		       double vB[3][3],double xB[3],double dxB[3])
{
  int iflag;
  int i,j,k;
  int i1,i2,j1,j2;
  double r,r0,r1;
  double d1,d2;
  double eps=1e-12;
  double D[3];
  double c[3][3];
  //
  // D=distance between centers
  // C=scalar product of axes
  //
  for(i=0;i<3;i++) D[i]=xB[i]-xA[i];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	c[i][j]=0;
	for(k=0;k<3;k++)
	  c[i][j]=c[i][j]+vA[i][k]*vB[j][k];
      }
  //
  // separating axes based on the faces of box A
  //
  for(i=0;i<3;i++)
    {
      r0=dxA[i];
      r1=0;
      r=0;
      for(j=0;j<3;j++)
	{
	  r1+=dxB[j]*fabs(c[i][j]);
	  r+=fabs(vA[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // separating axes based on the faces of box B
  //
  for(i=0;i<3;i++)
    {
      r1=dxB[i];
      r0=0;
      r=0;
      for(j=0;j<3;j++)
	{
	  r0+=dxA[j]*fabs(c[j][i]);
	  r+=fabs(vB[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // cross products
  //
  for(i=0;i<3;i++)
    {
      i1=(i+1)%3;
      i2=(i+2)%3;
      for(j=0;j<3;j++)
	{
	  j1=(j+1)%3;
	  j2=(j+2)%3;

	  r0=dxA[i1]*fabs(c[i2][j])+dxA[i2]*fabs(c[i1][j]);
	  r1=dxB[j1]*fabs(c[i][j2])+dxB[j2]*fabs(c[i][j1]);

	  d2=0;
	  d1=0;
	  for(k=0;k<3;k++)
	    {
	      d2+=vA[i2][k]*D[k];
	      d1+=vA[i1][k]*D[k];
	    }

	  r=fabs(c[i1][j]*d2-c[i2][j]*d1);

	  if (r > (r0+r1+eps)) {
	    return 0;
	  }
	}
    }
  //
  // return zero if no separation can be found
  //
  return 1;
}



void writebbox(OBB *obb,int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"qbox%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,
	  1);

  for(l=0;l<2;l++)
    {
      il=2*(l%2)-1;
      for(k=0;k<2;k++)
	{
	  ik=2*(k%2)-1;
	  for(j=0;j<2;j++)
	    {
	      ij=2*(j%2)-1;
	      xx[0]=xx[1]=xx[2]=0;
	      for(m=0;m<3;m++)
		xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
		  +ik*obb->vec[1][m]*obb->dxc[1]
		  +il*obb->vec[2][m]*obb->dxc[2];
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	    }
	}
    }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fclose(fp);
}



void writePoints(double *x,int nsearch,int bid)
{
  FILE *fp;
  char intstring[7];
  char fname[80];
  int i;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"points%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  for(i=0;i<nsearch;i++)
    fprintf(fp,"%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);
  fclose(fp);
}


/* ---- Linked-List Functions ---- */

void deallocateLinkList(DONORLIST *temp)
{
  if (temp!=NULL)
    {
      deallocateLinkList(temp->next);
      free(temp);
    }
}

void deallocateLinkList2(INTEGERLIST *temp)
{
  if (temp!=NULL)
    {
      deallocateLinkList2(temp->next);
      free(temp);
    }
}


void insertInList(DONORLIST **donorList,DONORLIST *temp1)
{
  DONORLIST *temp;
  DONORLIST *ptemp;
  int inserted;
  temp=*donorList;
  inserted=0;
  ptemp=NULL;
  while(temp!=NULL && !inserted)
    {
      if (fabs(temp->donorRes) > temp1->donorRes)
  {
    temp1->next=temp;
    if (ptemp!=NULL)
      {
        ptemp->next=temp1;
      }
    else
      {
        *donorList=temp1;
      }
    inserted=1;
  }
      else
  {
    ptemp=temp;
    temp=temp->next;
  }
    }
  if (!inserted)
    {
     if (*donorList)
      {
       temp1->next=NULL;
       ptemp->next=temp1;
      }
     else
       {
        temp1->next=NULL;
        *donorList=temp1;
       }
    }
}

void solvec(double **a,double *b,int *iflag,int n)
{
  int i,j,k,l,flag,temp1;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;


  for(i=0;i<n;i++)
  {
    if (fabs(a[i][i]) < eps)
    {
      flag=1;
      for(k=i+1;k<n && flag;k++)
      {
        if (a[k][i]!=0)
        {
          flag=0;
          for(l=0;l<n;l++)
          {
            temp=a[k][l];
            a[k][l]=a[i][l];
            a[i][l]=temp;
          }
          temp=b[k];
          b[k]=b[i];
          b[i]=temp;
        }
      }
      if (flag) {*iflag=0;return;}
    }
    for(k=i+1;k<n;k++)
    {
      if (i!=k)
      {
        fact=-a[k][i]/a[i][i];
        for(j=0;j<n;j++)
        {
          a[k][j]+=fact*a[i][j];
        }
        b[k]+=fact*b[i];
      }
    }
  }

  for(i=n-1;i>=0;i--)
  {
    sum=0;
    for(j=i+1;j<n;j++)
      sum+=a[i][j]*b[j];
    b[i]=(b[i]-sum)/a[i][i];
  }
  *iflag=1;
  return;

}

void newtonSolve(double f[7][3],double *u1,double *v1,double *w1)
{
  int i,j,k;
  int iter,itmax,isolflag;
  double u,v,w;
  double uv,wu,vw,uvw,norm,convergenceLimit;
  double *rhs;
  double **lhs;
  double alph;
  //
  lhs=(double **)malloc(sizeof(double)*3);
  for(i=0;i<3;i++)
    lhs[i]=(double *)malloc(sizeof(double)*3);
  rhs=(double *)malloc(sizeof(double)*3);
  //
  itmax=500;
  convergenceLimit=1e-14;
  alph=1.0;
  isolflag=1.0;
  //
  u=v=w=0.5;
  //
  for(iter=0;iter<itmax;iter++)
    {
      uv=u*v;
      vw=v*w;
      wu=w*u;
      uvw=u*v*w;

      for(j=0;j<3;j++)
  rhs[j]=f[0][j]+f[1][j]*u+f[2][j]*v+f[3][j]*w+
    f[4][j]*uv + f[5][j]*vw + f[6][j]*wu +
    f[7][j]*uvw;

      norm=rhs[0]*rhs[0]+rhs[1]*rhs[1]+rhs[2]*rhs[2];
      if (sqrt(norm) <= convergenceLimit) break;

      for(j=0;j<3;j++)
  {
    lhs[j][0]=f[1][j]+f[4][j]*v+f[6][j]*w+f[7][j]*vw;
    lhs[j][1]=f[2][j]+f[5][j]*w+f[4][j]*u+f[7][j]*wu;
    lhs[j][2]=f[3][j]+f[6][j]*u+f[5][j]*v+f[7][j]*uv;
  }

      solvec(lhs,rhs,&isolflag,3);
      if (isolflag==0) break;

      u-=(rhs[0]*alph);
      v-=(rhs[1]*alph);
      w-=(rhs[2]*alph);
    }

  if (isolflag==0) {
    u=1.0;
    v=w=0.;
  }
  *u1=u;
  *v1=v;
  *w1=w;
  for(i=0;i<3;i++) free(lhs[i]);
  free(lhs);
  free(rhs);
  return;
}


void computeNodalWeights(double xv[8][3],double *xp,double frac[8],int nvert)
{
  int i,j,k,isolflag;
  double **lhs;
  double *rhs;
  double f[8][3];
  double u,v,w;
  double oneminusU,oneminusV,oneminusW,oneminusUV;

  switch(nvert)
  {
    case 4:
      //
      // tetrahedron
      //
      lhs=(double **)malloc(sizeof(double)*3);
      for(i=0;i<3;i++)
        lhs[i]=(double *)malloc(sizeof(double)*3);
      rhs=(double *)malloc(sizeof(double)*3);
      for(k=0;k<3;k++)
      {
        for(j=0;j<3;j++)
          lhs[j][k]=xv[k][j]-xv[3][j];
        rhs[k]=xp[k]-xv[3][k];
      }
      //
      // invert the 3x3 matrix
      //
      solvec(lhs,rhs,&isolflag,3);
      //
      // check if the solution is not degenerate
      //
      if (isolflag)
      {
        for(k=0;k<3;k++) frac[k]=rhs[k];
        frac[3]=1.-frac[0]-frac[1]-frac[2];
      }
      else
      {
        frac[0]=1.0;
        frac[1]=frac[2]=frac[3]=0;
      }
      for(i=0;i<3;i++) free(lhs[i]);
      free(lhs);
      free(rhs);
      break;

    case 5:
      //
      // pyramid
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[3][j]-xv[0][j];
        f[3][j]=xv[4][j]-xv[0][j];

        f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
        f[5][j]=xv[0][j]-xv[3][j];
        f[6][j]=xv[0][j]-xv[1][j];
        f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j];
      }

      newtonSolve(f,&u,&v,&w);
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;

      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=w;

      break;

    case 6:
      //
      // prizm
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[2][j]-xv[0][j];
        f[3][j]=xv[3][j]-xv[0][j];
        //
        f[4][j]=0;
        f[5][j]=xv[0][j]-xv[2][j]-xv[3][j]+xv[5][j];
        f[6][j]=xv[0][j]-xv[1][j]-xv[3][j]+xv[4][j];
        f[7][j]=0.;
      }

      newtonSolve(f,&u,&v,&w);

      oneminusUV=1.0-u-v;
      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;

      frac[0]=oneminusUV*oneminusW;
      frac[1]=u*oneminusW;
      frac[2]=v*oneminusW;
      frac[3]=oneminusUV*w;
      frac[4]=u*w;
      frac[5]=v*w;

      break;

    case 8:
      //
      // hexahedra
      //
      for(j=0;j<3;j++)
      {
        f[0][j]=xv[0][j]-xp[j];
        f[1][j]=xv[1][j]-xv[0][j];
        f[2][j]=xv[3][j]-xv[0][j];
        f[3][j]=xv[4][j]-xv[0][j];

        f[4][j]=xv[0][j]-xv[1][j]+xv[2][j]-xv[3][j];
        f[5][j]=xv[0][j]-xv[3][j]+xv[7][j]-xv[4][j];
        f[6][j]=xv[0][j]-xv[1][j]+xv[5][j]-xv[4][j];
        f[7][j]=-xv[0][j]+xv[1][j]-xv[2][j]+xv[3][j]+
            xv[4][j]-xv[5][j]+xv[6][j]-xv[7][j];
      }

      newtonSolve(f,&u,&v,&w);

      oneminusU=1.0-u;
      oneminusV=1.0-v;
      oneminusW=1.0-w;

      frac[0]=oneminusU*oneminusV*oneminusW;
      frac[1]=u*oneminusV*oneminusW;
      frac[2]=u*v*oneminusW;
      frac[3]=oneminusU*v*oneminusW;
      frac[4]=oneminusU*oneminusV*w;
      frac[5]=u*oneminusV*w;
      frac[6]=u*v*w;
      frac[7]=oneminusU*v*w;

      break;

    default:
      printf("Interpolation not implemented for polyhedra with %d vertices\n",nvert);
      break;
  }
}

double computeCellVolume(double xv[8][3],int nvert)
{
 double vol;
 int itype;
 int nfaces;
 int numverts[4][6]={3,3,3,3,0,0,4,3,3,3,3,0,3,4,4,4,3,0,4,4,4,4,4,4};
 int faceInfo[4][24]={1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0,
                       1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0,
                       1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0,
                       1,2,3,4,1,5,6,2,2,6,7,3,3,7,8,4,1,4,8,5,5,8,7,6};
 switch(nvert)
 {
   case 4:
     itype=0;
     nfaces=4;
     break;
   case 5:
     itype=1;
     nfaces=5;
     break;
   case 6:
     itype=2;
     nfaces=5;
     break;
   case 8:
     itype=3;
     nfaces=6;
     break;
   default:
     printf("Invalid value of nvert (%d) in computeCellVolume!\n",nvert);
     exit(1);
 }

 cellvolume_(&vol,xv,&numverts[itype],&faceInfo[itype],&nfaces,&nvert);
 return vol;
}
