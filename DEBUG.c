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
	printArrayDouble(KM->cluster_radius, KM->k, "radius of cluster");

	printf("cluster_centroid: \n");
	printArraysDouble(KM->cluster_centroid, KM->k, KM->dim, "centroid of cluster");

	printf("cluster_assign: \n");
	printArraysInt(KM->cluster_assign, KM->ndata, "cluster of DP");

	// printf("cluster_group: \n");
	// printArraysInt(KM->cluster_group, KM->ndata, "data point-> ");

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

void writeResults(int dim, int ndata, double* data, int* cluster_assign)
{
	int i;
	FILE* file = fopen("data.txt", "w");

	fprintf(file, "%d\n", dim);
	fprintf(file, "%d\n", ndata);
	for (i = 0; i < dim * ndata; i++)
		fprintf(file, "%lf\n", data[i]);

	for (i = 0; i < ndata; i++)
		fprintf(file, "%d\n", cluster_assign[i]);

	fclose(file);
}

void printArray(int * nums, int count)
{
	int i;
	for (i = 0; i < count; i++)
	{
		printf("%d ", nums[i]);
	}
	printf("\n");
}

void printArrayDoubles(double * nums, int ndata, int dim)
{
	int i,j,first_index;
	for (i = 0; i < ndata; i++)
	{
		first_index = i * dim;
		for (j = 0; j < dim; j++)
		{
			printf("%lf ", nums[first_index + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printArrayDouble(double * ArrayDouble, int size, const char * text)
{
	int i, j;
	for (i = 0; i < size; i++)
	{
		printf("%s %d = ( ", text, i);
		for (j = 0; j < 1; j++)
		{
			if (j != 0)
			{
				printf("%lf,", ArrayDouble[i]);
			}
			else
			{
				printf("%lf", ArrayDouble[i]);
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

void printStack(struct stackBase *stack)
{
	struct stackNode * iterator = stack->firstNode;
	int i;
	printf("Stack: \n");
	printf("stackDepth = %d, firstNode = %p \n", stack->stackDepth,(void *)stack->firstNode);
	while (iterator != NULL)
	{
		for (i = 0; i < stack->arraySize; i++)
		{
			printf("dim: %d = %lf, \n", i, (iterator->pointArray)[i]);
		}
		printf("distance = %lf, cluster = %d ", iterator->distance, iterator->cluster);
		printf("\n\n");
		iterator = iterator->nextNode;
	}
}
