#include "mainHeader.h"


void initializeKM(struct kmeans ** KM, int dim, int ndata, int subdomain, double * dataArray, int k)
{
  *KM = (struct kmeans *)malloc(sizeof(struct kmeans));
  (*KM)->dim = dim;
  (*KM)->ndata = ndata;
  (*KM)->subdomain = subdomain;
  (*KM)->data = dataArray;
  (*KM)->k = k;
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
		target[i] = 0;
	}
	return target;
}

void kmeans(struct kmeans ** KM, int dim, int ndata, int subdomain, double * dataArray, int k)
{
  initializeKM(KM,dim,ndata,subdomain,dataArray,k);
  GetKCentroids(*KM);
  ClusterizeKM(*KM);
  return;
}
