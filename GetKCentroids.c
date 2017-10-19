#include "mainHeader.h"

void GetKCentroids(struct kmeans * KM)
{
  /* 1. Randomly select a point in the data array for the first centroid. */

  int i,size = KM->dim;
  double * first_cent = (double *)malloc(sizeof(double) * size);
  if (KM->world_rank == 0)
  {
    int first_centroid = 0;
    for (i = 0; i < size; i++)
    {
      first_cent[i] = (KM->data)[(first_centroid * (size)) + i];
    }
  }
  MPI_Bcast(first_cent,size,MPI_DOUBLE,0,MCW);
  MPI_Barrier(MCW);
  for (i = 0; i < size; i++)
  {
    (KM->cluster_centroid)[0][i] = first_cent[i];
  }

  /* 2. Select each subsequent centroid based on the max-of-the-min criteria given in class. */

  int numClusters = 1;
  while (numClusters < KM->k)
  {
    numClusters = GetNextCluster(KM,numClusters);
  }
  return;
}

int GetNextCluster(struct kmeans * KM, int numCentroids)
{
  double * nextCent = (double *)malloc(sizeof(double) * (KM->dim + 1)); //First index stores maxminDist of point.
  double * AllNextCentroids = (double *)malloc(sizeof(double) * ((KM->dim + 1) * KM->world_size));
  int i,j,nextCentroid,first_index;
  double * minDistArray = allocateAndInitializeZeroDouble((KM->ndata));
  double minDist = INFINITY,maxminDist = -INFINITY, distance;

  for (i = 0; i < KM->ndata; i++)
  {
    first_index = (KM->dim) * i;
    minDist = INFINITY;
    for (j = 0; j < numCentroids; j++)
    {
      distance = GetDistance2PointsDC(KM,first_index,j);
      if (distance < minDist)
      {
        minDist = distance;
        minDistArray[i] = distance;
      }
    }
  }

  for (i = 0; i < KM->ndata; i++)
  {
    if (minDistArray[i] > maxminDist)
    {
      maxminDist = minDistArray[i];
      nextCentroid = i;
    }
  }
  for (i = 1; i < (KM->dim + 1); i++)
  {
    nextCent[i] = (KM->data)[(nextCentroid * KM->dim) + (i-1)];
  }
  nextCent[0] = maxminDist;

  MPI_Allgather(nextCent,(KM->dim + 1),MPI_DOUBLE,AllNextCentroids,(KM->dim + 1),MPI_DOUBLE,MCW);
  MPI_Barrier(MCW);

  double maxWorldDist = -INFINITY;
  int stride = 0,minLoc;
  for (i = 0; i < KM->world_size; i++)
  {
    stride = i * (KM->dim + 1);
    if (AllNextCentroids[stride] > maxWorldDist)
    {
      maxWorldDist = AllNextCentroids[stride];
      minLoc = stride+1;
    }
  }

  for (i = 0; i < KM->dim; i++)
  {
    (KM->cluster_centroid)[numCentroids][i] = AllNextCentroids[minLoc + i];
  }

  free(minDistArray);
  free(nextCent);
  free(AllNextCentroids);
  return numCentroids + 1;
}

double GetDistance2PointsDC(struct kmeans *KM, int first_index, int centroid)
{
  int i;
  double distCalc = 0;
  for (i = 0; i < KM->dim; i++)
  {
    distCalc += pow( fabs((KM->data)[first_index + i] - (KM->cluster_centroid)[centroid][i]), 2);
  }
  distCalc = sqrt(distCalc);
  return distCalc;
}
