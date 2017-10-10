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
  if (FIRST_CENTROID) { printf("world_rank %d\n",KM->world_rank); printArrayDouble(first_cent,size,"value->"); }
  for (i = 0; i < size; i++)
  {
    (KM->cluster_centroid)[0][i] = first_cent[i];
  }
  if (FIRST_CENTROID) { printf("world_rank %d\n",KM->world_rank); printArrayDouble((KM->cluster_centroid)[0],size,"value->"); }

  /* 2. Select each subsequent centroid based on the max-of-the-min criteria given in class. */

  int numClusters = 1;
  if (DEBUG_SELECTK) { printf("finished\n");
                        printf("centroids:\n");
                        printArraysDouble(KM->cluster_centroid,KM->k,KM->dim,"centroid -> "); }
  while (numClusters < KM->k)
  {
    numClusters = GetNextCluster(KM,numClusters);
    if (WAYPOINTS) { printf("%d centroids obtained.\n",numClusters);}
  }
  if (NEXT_CENTROID1) { printf("finished\n");
                        printf("centroids in world_rank %d:\n",KM->world_rank);
                        printArraysDouble(KM->cluster_centroid,KM->k,KM->dim,"centroid -> "); }
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
      if (DEBUG_SELECTK) printf("distance between %d and centroid %d = %lf \n",i,j,distance);
      if (distance < minDist)
      {
        minDist = distance;
        minDistArray[i] = distance;
      }
    }
  }

  if (DEBUG_SELECTK) { printf("minDistArray:\n"); printArrayKMD(minDistArray,KM->ndata); }
  for (i = 0; i < KM->ndata; i++)
  {
    if (minDistArray[i] > maxminDist)
    {
      maxminDist = minDistArray[i];
      nextCentroid = i;
    }
  }
  if (DEBUG_SELECTK) printf("minDistArray[%d] = %lf, maxminDist = %lf \n",nextCentroid,minDistArray[nextCentroid],maxminDist);
  for (i = 1; i < (KM->dim + 1); i++)
  {
    nextCent[i] = (KM->data)[(nextCentroid * KM->dim) + (i-1)];
  }
  nextCent[0] = maxminDist;

  if (NEXT_CENTROID) { printf("world_rank %d\n",KM->world_rank); printArrayDouble(nextCent,KM->dim + 1,"nextCent->"); }
  MPI_Allgather(nextCent,(KM->dim + 1),MPI_DOUBLE,AllNextCentroids,(KM->dim + 1),MPI_DOUBLE,MCW);
  MPI_Barrier(MCW);
  if (NEXT_CENTROID) { printf("world_rank %d\n",KM->world_rank); printArrayDouble(AllNextCentroids,((KM->dim + 1) * KM->world_size),"value->"); }

  double minWorldDist = INFINITY;
  int stride = 0,minLoc;
  for (i = 0; i < KM->world_size; i++)
  {
    stride = i * (KM->dim + 1);
    if (AllNextCentroids[stride] < minWorldDist)
    {
      minWorldDist = AllNextCentroids[stride];
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
