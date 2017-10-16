#include "mainHeader.h"

void GetKCentroids(struct kmeans * KM)
{
  /* 1. Randomly select a point in the data array for the first centroid. */
  int first_centroid = 0;
  int i;
  for (i = 0; i < KM->dim; i++)
  {
    (KM->cluster_centroid)[0][i] = (KM->data)[(first_centroid * (KM->dim)) + i];
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
  for (i = 0; i < KM->dim; i++)
  {
    (KM->cluster_centroid)[numCentroids][i] = (KM->data)[(nextCentroid * KM->dim) + i];
  }
  free(minDistArray);
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
