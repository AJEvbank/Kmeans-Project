#include "mainHeader.h"

void quickSort( struct kmeans * KM, int l, int r)
{

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

int partition( struct kmeans * KM, int left, int right)
{

  int pivot, i, j, temp, jump;
  pivot = (KM->cluster_assign)[left];
  i = left; j = right + 1;

  while( 1)
  {
  	do ++i; while( (KM->cluster_assign)[i] <= pivot && i <= right );
  	do --j; while( (KM->cluster_assign)[j] > pivot );
  	if( i >= j ) { jump = j; break; }
    /* Swap data points here. */
    SwapPoints(KM,i,j);
  	temp = (KM->cluster_assign)[i];
    (KM->cluster_assign)[i] = (KM->cluster_assign)[j];
    (KM->cluster_assign)[j] = temp;
  }
  /* Swap the pivot point. */
  SwapPoints(KM,left,jump);
  temp = (KM->cluster_assign)[left];
  (KM->cluster_assign)[left] = (KM->cluster_assign)[jump];
  (KM->cluster_assign)[jump] = temp;
  return jump;
}

void SwapPoints(struct kmeans * KM,int left,int right)
{
  double temp;
  int i,left_index = left * (KM->dim),right_index = right * (KM->dim);
  for (i = 0; i < (KM->dim); i++)
  {
    temp = (KM->data)[left_index + i];
    (KM->data)[left_index + i] = (KM->data)[right_index + i];
    (KM->data)[right_index + i] = temp;
  }
  return;
}
