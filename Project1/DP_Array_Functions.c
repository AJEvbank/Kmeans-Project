#include "myHeader.h"

double * dataArrayGenerator(int dim, int ndata, double max_double)
{
	int i = 0, first_dim_of_data_point;

	double * dataArray = (double *)malloc(dim * ndata * sizeof(double));

	for (i = 0; i < ndata; i++)
	{
		first_dim_of_data_point = dim * i;
		int j = 0;
		for (j = 0; j < dim; j++)
		{
			dataArray[first_dim_of_data_point + j] = (((double)rand()/(double)(RAND_MAX)) * (max_double));			
		}			
	}
	return dataArray;
}

int bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * result)
{
	int first_index, i, j, nearestPoint,isQuery = 1;
	double minDist = INFINITY, calcDist = 0;
	for (i = 0; i < ndata; i++)
	{
		calcDist = 0;
		first_index = i * dim;
		for (j = 0; j < dim; j++)
		{
			calcDist += pow(fabs(query[j] - dataArray[first_index + j]), 2);
		}
		calcDist = sqrt(calcDist);
		if (calcDist < minDist)
		{
			minDist = calcDist;
			nearestPoint = i;
		}
	}
	printf("Brute force computation: \n");
	printf("Data Point: %d/%d, distance = %lf \n", nearestPoint, ndata, minDist);
	printf("[ ");
	for (i = 0; i < dim; i++)
	{
		printf("%lf ", dataArray[(nearestPoint*dim) + i]);
		if (dataArray[(nearestPoint*dim) + i] != result[i])
		{
			isQuery = 0;
		}
	}
	printf("]\n");
	
	return isQuery;
}


