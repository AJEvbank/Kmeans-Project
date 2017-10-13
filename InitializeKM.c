#include "mainHeader.h"


void initializeKM(struct kmeans ** KM, int dim, int ndata, double * dataArray, int k, int world_rank, int world_size)
{
  *KM = (struct kmeans *)malloc(sizeof(struct kmeans));
  (*KM)->dim = dim;
  (*KM)->ndata = ndata;
  (*KM)->data = dataArray;
  (*KM)->k = k;
  (*KM)->world_rank = world_rank;
  (*KM)->world_size = world_size;
  (*KM)->cluster_size = allocateAndInitializeZeroInt(k);
  (*KM)->cluster_start = allocateAndInitializeZeroInt(k);
  (*KM)->cluster_radius = allocateAndInitializeZeroDouble(k);
  (*KM)->cluster_assign = allocateAndInitializeZeroInt(ndata);
  (*KM)->cluster_centroid = allocateAndInitializeZeroDoubleMulti(k,dim);
  //(*KM)->cluster_group = allocateAndInitializeZeroInt(ndata);
  return;
}

int * allocateAndInitializeZeroInt(int size_of_target)
{
	int * target = (int *)malloc(size_of_target * sizeof(int));
	int i = 0;
	for (i = 0; i < size_of_target; i++) {
		target[i] = 0;
	}
	return target;
}

double ** allocateAndInitializeZeroDoubleMulti(int k, int dimension)
{
	double ** target = (double **)malloc(sizeof(double *)* k);
	int i,j;
	for (i = 0; i < k; i++)
	{
		target[i] = (double *)malloc(sizeof(double) * dimension);
		for (j = 0; j < dimension; j++)
		{
			target[i][j] = 0;
		}
	}
	return target;
}

double * allocateAndInitializeZeroDouble(int size_of_target)
{
	double * target = (double *)malloc(size_of_target * sizeof(double));
	int i = 0;
	for (i = 0; i < size_of_target; i++) {
		target[i] = 0.0;
	}
	return target;
}

void kmeans(struct kmeans ** KM, int dim, int ndata, double * dataArray, int k, int world_rank, int world_size)
{
  initializeKM(KM,dim,ndata,dataArray,k,world_rank,world_size);
  if (WAYPOINTS2) { printf("Kmeans initialized on world_rank %d.\n",(*KM)->world_rank); }
  GetKCentroids(*KM);
  MPI_Barrier(MCW);
  if (WAYPOINTS) { printf("Got k centroids on world_rank %d.\n",(*KM)->world_rank); }
  ClusterizeKM(*KM, THRESHOLD);
  MPI_Barrier(MCW);
  if (WAYPOINTS2) printf("Clusterized on world_rank %d \n",(*KM)->world_rank);
  return;
}
