#ifndef _CODETYPES_H
#define _CODETYPES_H

#define MPICH_SKIP_MPICXX
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include "math.h"
#include "mpi.h"

using std::min;
using std::max;
using std::abs;

/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double

/*====================================================================*/
/*  Base for indexing (0 or 1)
/*====================================================================*/
# define BASE 0

/*====================================================================*/
/*  Define TIOGA conventions for node/cell blanking
/*====================================================================*/
#define NORMAL 1
#define HOLE 0
#define FRINGE -1

/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
#define ZERO               0.0e+00
#define ONE                1.0e+00
#define TWO                2.0e+00
#define THREE              3.0e+00
#define FOUR               4.0e+00
#define HALF               0.5e+00
#define THIRD              0.333333333e+00
#define FIFTH              0.2
#define PI                 3.1415926535897932e+00
#define RAD2DEG            (180.0/PI)
#define DEG2RAD            (PI/180.0)
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define NFRINGE            3
#define NVAR               6

/*==================================================================*/
/* inline debugging tools                                             */
/*==================================================================*/
# define tracei(x)  printf("#tioga:\t"#x" =%d\n",x);
# define traced(x)  printf("#tioga:\t"#x" =%.16e\n",x);
//# define min(x,y)  (x) < (y) ? (x) : (y)
//# define max(x,y)  (x) > (y) ? (x) : (y)
# define debug(x,y)  printf("#tioga:\t"#x"=%d,"#y"=%d\n",x,y);
# define stdwrite(x) if (myid==0) printf("#tioga:\t"#x"\n");
# define dstr(x) printf("#tioga:\t"#x"\n");
# define ditch(x,y) {dstr(x); tracei(y); MPI_Abort(MPI_COMM_WORLD,ierr);}

/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
#define Sign(a1,a2)\
        (((a2) < ZERO)? - fabs(a1): fabs(a1))
#define Max(a1,a2)\
        (((a1) >= (a2))? (a1): (a2))
#define Min(a1,a2)\
        (((a1) <= (a2))? (a1): (a2))
#define Abs(aa)\
        (((aa) >= 0)? (aa): -(aa))
#define Round(aa)\
        (int) ((fabs((aa) - floor(aa)) >= HALF)? ceil(aa): floor(aa))
//#define swap(a,b) { a=a+b;b=a-b;a=a-b;}

/*********************************************************************/
/* Code specific types
/*********************************************************************/
typedef struct HOLEMAP
{
  int existWall;
  int nx[3];
  int *samLocal;
  int *sam;
  double extents[6];
} HOLEMAP;


typedef struct OBB
{
  double xc[3];
  double dxc[3];
  double vec[3][3];
} OBB;

typedef struct DONORLIST
{
  int donorData[3];
  double donorRes;
  struct DONORLIST *next;
} DONORLIST;

typedef struct PACKET
{
  int nints;
  int nreals;
  int *intData;
  REAL *realData;
} PACKET;

typedef struct INTERPLIST
{
  int cancel;
  int nweights;
  int receptorInfo[2];
  int *inode;
  double *weights;
} INTERPLIST;

typedef struct INTEGERLIST
{
  int inode;
  struct INTEGERLIST *next;
} INTEGERLIST;

#endif // _CODETYPES_H
