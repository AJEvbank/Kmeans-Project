#include "mainHeader.h"

void displayKM(struct kmeans * KM)
{
	printf("Kmeans: \n");
	printf("dim = %d, ndata = %d, k = %d \n", KM->dim, KM->ndata, KM->k);

	printf("Working data array: \n");
	printDataArray(KM->data, KM->dim, KM->ndata);

	printf("cluster_size: \n");
	printArraysInt(KM->cluster_size, KM->k, "size of cluster");

	printf("cluster_start: \n");
	printArraysInt(KM->cluster_start, KM->k, "start of cluster");

	printf("cluster_radius: \n");
	printArrayKMD(KM->cluster_radius, KM->k);

	printf("cluster_centroid: \n");
	printArraysDouble(KM->cluster_centroid, KM->k, KM->dim, "centroid of cluster");

	printf("cluster_assign: \n");
	printArraysInt(KM->cluster_assign, KM->ndata, "cluster of DP");

	return;
}

void printArraysInt(int * ArrayInt, int size, const char * text)
{
	int i;
	for (i = 0; i < size; i++)
	{
		printf("%s %d = %d, \n", text, i, ArrayInt[i]);
	}
	return;
}

void printArraysDouble(double ** ArrayDouble, int size, int dim, const char * text)
{
	int i, j;
	for (i = 0; i < size; i++)
	{
		printf("%s %d = ( ", text, i);
		for (j = 0; j < dim; j++)
		{
			if (j != dim - 1)
			{
				printf("%lf,", ArrayDouble[i][j]);
			}
			else
			{
				printf("%lf", ArrayDouble[i][j]);
			}
		}
		if (i != size - 1)
		{
			printf("),\n");
		}
		else
		{
			printf(")\n");
		}
	}
	return;
}

void printArrayKMD(double * array, int size)
{
	int i;
	printf("[ ");
	for (i = 0; i < size; i++)
	{
		printf("%lf ", array[i]);
	}
	printf("]\n");
}

void printDataArray(double * dataArray, int dim, int ndata)
{
	int i, j;
	for (i = 0; i < ndata; i++)
	{
		int first_dim_of_data_point = dim * i;
		for (j = 0; j < dim; j++)
		{
			printf("%lf, ", dataArray[first_dim_of_data_point + j]);

		}
		printf("\n");
	}
	return;
}
