#include "codetypes.h"

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

