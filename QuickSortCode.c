#include "mainHeader.h"

void quickSort( struct kmeans * KM, int l, int r)
{

  int * a = KM->cluster_assign;
  int j;

  if( l < r )
  {
  	// divide and conquer
      j = partition( KM, l, r);
     quickSort( KM, l, j-1);
     quickSort( KM, j+1, r);
  }
  return;
}



int partition( struct kmeans * KM, int l, int r)
{

  int * a = KM->cluster_assign;
  int pivot, i, j, t;
  pivot = a[l];
  i = l; j = r+1;

  while( 1)
  {
  	do ++i; while( a[i] <= pivot && i <= r );
  	do --j; while( a[j] > pivot );
  	if( i >= j ) break;
    /* Swap data points here. */
    SwapPoints(KM,i,j);
  	t = a[i]; a[i] = a[j]; a[j] = t;
  }
  SwapPoints(KM,l,j);
  t = a[l]; a[l] = a[j]; a[j] = t;
  return j;
}

void SwapPoints(KM,int left,int right)
{
  double temp;
  int i;
  for (i = 0; i < KM->dim; i++)
  {
    temp = (KM->data)[left + i];
    (KM->data)[left + i] = (KM->data)[right + i];
    (KM->data)[right + i] = temp;
  }
  return;
}
