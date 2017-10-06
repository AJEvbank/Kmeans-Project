#include "mainHeader.h"

void ClusterizeKM(struct kmeans * KM, int threshold)
{
  if (DEBUG_THRESHOLD)
  {
    printf("Centroids before first iteration \n");
    printArraysDouble(KM->cluster_centroid,KM->k,KM->dim,"centroid");
    printf("\n\n");
  }
  int i,changed;
  for (i = 0,changed = 1; changed != 0; i++)
  {

  /* Assign each data point to the cluster with the nearest centroid */
  AssignDPs(KM);

  /* For each cluster, recalculate the centroid based on the data oints newly assigned. */
  changed = RecalculateCentroids(KM);

  /* Repeat the iteration until cluster assignments do not change or threshold is reached. */

    if (DEBUG_THRESHOLD)
    {
      printf("Centroids on iteration %d/%d \n",i,threshold);
      printArraysDouble(KM->cluster_centroid,KM->k,KM->dim,"centroid");
      if (changed)
      {
        printf("Changes occurred...");
      }
      else
      {
        printf("No changes occurred...");
      }
      printf("\n\n");
    }
  }
  SaveClusters(KM);

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

int RecalculateCentroids(struct kmeans * KM)
{
  int * newClusterSizes = allocateAndInitializeZeroInt(KM->k);
  double * centroids = allocateAndInitializeZeroDouble(KM->k * KM->dim);
  int i,j,group,group_coordinate,first_index,changed = 0;

  /* Get the sum of each dimension in each cluster */
  for ( i = 0; i < KM->ndata; i++)
  {
    group = (KM->cluster_assign)[i];
    newClusterSizes[group] += 1;
    group_coordinate = group * KM->dim;
    first_index = i * KM->dim;
    for (j = 0; j < KM->dim; j++)
    {
      centroids[group_coordinate + j] += (KM->data)[first_index + j];
    }
  }

  /* Divide the sum of each dimension in each cluster by the size of the cluster. */
  for ( i = 0; i < KM->k; i++)
  {
    first_index = i * KM->dim;
    for ( j = 0; j < KM->dim; j++)
    {
      centroids[first_index + j] /= (double)newClusterSizes[i];
    }
  }

  /* Save the newly calculated centroids in the kmeans structure. */
  for ( i = 0; i < KM->k; i++)
  {
    first_index = i * KM->dim;
    for ( j = 0; j < KM->dim; j++)
    {
      if ((KM->cluster_centroid)[i][j] != centroids[first_index + j])
      {
        (KM->cluster_centroid)[i][j] = centroids[first_index + j];
        changed = 1;
      }
    }
  }

  /* Save the newly calculated cluster sizes in the kmeans structure. */
  for ( i = 0; i < KM->k; i++)
  {
    if ((KM->cluster_size)[i] != newClusterSizes[i])
    {
      (KM->cluster_size)[i] = newClusterSizes[i];
      changed = 1;
    }
  }

  free(newClusterSizes);
  free(centroids);
  return changed;
}

void SaveClusters(struct kmeans * KM)
{
  /* Sort the data array (locally) by cluster assignment */

  SortDataArray(KM);

  /* Calculate the radius and start of each cluster. */


  return;
}

void SortDataArray(struct kmeans * KM)
{


  return;
}
