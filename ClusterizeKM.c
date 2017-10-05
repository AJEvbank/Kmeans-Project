#include "mainHeader.h"

void ClusterizeKM(struct kmeans * KM, int threshold)
{
  int i;
  for (i = 0; i < threshold; i++)
  {

  /* Assign each data point to the cluster with the nearest centroid */
  AssignDPs(KM);

  /* For each cluster, recalculate the centroid based on the data oints newly assigned. */


  /* Repeat the iteration until cluster assignments do not change or threshold is reached. */
    return;
  }



  return;
}

void AssignDPs(struct kmeans * KM)
{
  int i,j,first_index,assignment;
  double distance = 0, minDist;
  for (i = 0; i < (KM->ndata); i++)
  {
    first_index = i * (KM->dim);
    minDist = INFINITY;
    for (j = 0; j < (KM->k); j++)
    {
      distance = GetDistance2PointsDC(KM, first_index, j);
      //if (DEBUG_ASSIGNA) printf("distance between %d and Centroid %d = %lf \n",i,j,distance);
      if (distance < minDist)
      {
        minDist = distance;
        assignment = j;
      }
    }
    (KM->cluster_assign)[i] = assignment;
  }

  return;
}
