#include "mainHeader.h"


void ClusterizeKM(struct kmeans * KM)
{

  int i,changed;
  for (i = 0,changed = 1; changed != 0; i++)
  {

  /* Assign each data point to the cluster with the nearest centroid */
  AssignDPs(KM);
  if (WAYPOINTS) { printf("Assigned at iteration %d.\n",i); }
  /* For each cluster, recalculate the centroid based on the data oints newly assigned. */
  changed = RecalculateCentroids(KM);
  if (WAYPOINTS) { printf("Centroids recalculated at iteration %d.\n",i); }

  /* Repeat the iteration until cluster assignments do not change or threshold is reached. */

  }
  SaveClusters(KM);
  if (WAYPOINTS) { printf("Clusters saved.\n"); }

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
  int * newClusterSizesLoc = allocateAndInitializeZeroInt(KM->k);
  int * newClusterSizesGlob = allocateAndInitializeZeroInt(KM->k);
  double * centroidsLoc = allocateAndInitializeZeroDouble(KM->k * KM->dim);
  double * centroidsGlob = allocateAndInitializeZeroDouble(KM->k * KM->dim);
  int i,j,group,group_coordinate,first_index,changed = 0;

  /* Get the sum of each dimension in each cluster */
  for ( i = 0; i < KM->ndata; i++)
  {
    group = (KM->cluster_assign)[i];
    newClusterSizesLoc[group] += 1;
    group_coordinate = group * KM->dim;
    first_index = i * KM->dim;
    for (j = 0; j < KM->dim; j++)
    {
      centroidsLoc[group_coordinate + j] += (KM->data)[first_index + j];
    }
  }

  /* Save the newly calculated cluster sizes in the kmeans structure. */
  for ( i = 0; i < KM->k; i++)
  {
    (KM->cluster_size)[i] = newClusterSizesLoc[i];
  }

  /* Get the global sums of each dimension and the global cluster sizes. */
  MPI_Allreduce(
                centroidsLoc,
                centroidsGlob,
                (KM->k * KM->dim),
                MPI_DOUBLE,
                MPI_SUM,
                MCW
  );
  MPI_Barrier(MCW);

  MPI_Allreduce(
                newClusterSizesLoc,
                newClusterSizesGlob,
                (KM->k),
                MPI_INT,
                MPI_SUM,
                MCW
  );
  MPI_Barrier(MCW);

  /* Divide the sum of each dimension in each cluster by the size of the cluster. */
  for ( i = 0; i < KM->k; i++)
  {
    first_index = i * KM->dim;
    for ( j = 0; j < KM->dim; j++)
    {
      centroidsGlob[first_index + j] /= (double)newClusterSizesGlob[i];
    }
  }

  /* Save the newly calculated centroids in the kmeans structure. */
  for ( i = 0; i < KM->k; i++)
  {
    first_index = i * KM->dim;
    for ( j = 0; j < KM->dim; j++)
    {
      if ((KM->cluster_centroid)[i][j] != centroidsGlob[first_index + j])
      {
        (KM->cluster_centroid)[i][j] = centroidsGlob[first_index + j];
        changed = 1;
      }
    }
  }

  free(newClusterSizesLoc);
  free(newClusterSizesGlob);
  free(centroidsLoc);
  free(centroidsGlob);
  return changed;
}

void SaveClusters(struct kmeans * KM)
{
  /* Sort the data array (locally) by cluster assignment */

  SortDataArray(KM);

  /* Set the cluster start values. */

  int i,count = 0;
  for ( i = 0; i < (KM->k); i++)
  {
    (KM->cluster_start)[i] = count;
    count += (KM->cluster_size)[i];
  }

  /* Calculate the radius of each cluster. */

  CalculateRadii(KM);

  /* Deal with empty clusters. */

  EmptyClusters(KM);

  return;
}

void SortDataArray(struct kmeans * KM)
{
  quickSort(KM, 0, (KM->ndata)-1);
  return;
}

void CalculateRadii(struct kmeans * KM)
{
  double * MaxRadiiLoc = allocateAndInitializeZeroDouble(KM->k);
  double * MaxRadiiGlob = allocateAndInitializeZeroDouble(KM->k);
  int i,j,start;
  double distance = 0;
  for (i = 0; i < (KM->k); i++)
  {
    for (j = 0; j < (KM->cluster_size)[i]; j++)
    {
      start = (KM->cluster_start)[i];
      distance = GetDistance2PointsDC(KM,(start + j)*KM->dim,i);
      if (distance > MaxRadiiLoc[i])
      {
        MaxRadiiLoc[i] = distance;
      }
    }
  }

  MPI_Allreduce(
                MaxRadiiLoc,
                MaxRadiiGlob,
                (KM->k),
                MPI_DOUBLE,
                MPI_MAX,
                MCW
  );
  MPI_Barrier(MCW);


  for (i = 0; i < (KM->k); i++)
  {
    (KM->cluster_radius)[i] = MaxRadiiGlob[i];
  }
  free(MaxRadiiLoc);
  free(MaxRadiiGlob);
  return;
}

void EmptyClusters(struct kmeans * KM)
{
  int * ClusterSizeCheckLoc = (int *)malloc(sizeof(int) * KM->k);
  int * ClusterSizeCheckGlob = (int *)malloc(sizeof(int) * KM->k);
  int i,limit = KM->k;
  for (i = 0; i < KM->k; i++)
  {
    ClusterSizeCheckLoc[i] = (KM->cluster_size)[i];
  }

  MPI_Allreduce(
                ClusterSizeCheckLoc,
                ClusterSizeCheckGlob,
                (KM->k),
                MPI_INT,
                MPI_SUM,
                MCW
  );
  MPI_Barrier(MCW);


  for (i = limit-1; i >= 0; i--)
  {
    if(ClusterSizeCheckGlob[i] == 0)
    {
      DeleteEmptyCluster(KM,i);
    }
  }
  free(ClusterSizeCheckLoc);
  free(ClusterSizeCheckGlob);
  return;
}

void DeleteEmptyCluster(struct kmeans * KM, int cluster)
{
  int i,limit = KM->k;

  /* Free the centroid. */
  free((KM->cluster_centroid)[cluster]);
  (KM->cluster_centroid)[cluster] = NULL;

  /* Shift back size, radius, start, and centroid. */
  for (i = cluster; i < limit-1; i++)
  {
    (KM->cluster_size)[i] = (KM->cluster_size)[i+1];
    (KM->cluster_radius)[i] = (KM->cluster_radius)[i+1];
    (KM->cluster_start)[i] = (KM->cluster_start)[i+1];
    (KM->cluster_centroid)[i] = (KM->cluster_centroid)[i+1];
  }

  /* Adjust the cluster assignment indices. */
  for (i = (KM->ndata); i >= (KM->cluster_start)[cluster-1]; i--)
  {
    if ((KM->cluster_assign)[i] > cluster)
    {
      (KM->cluster_assign)[i] = (KM->cluster_assign)[i] - 1;
    }
  }

  /* Decrement k. */
  KM->k = KM->k - 1;

  return;
}
