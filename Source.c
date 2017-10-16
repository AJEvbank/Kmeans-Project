#include "mainHeader.h"

int main(int argc, char** argv) {

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//MPI_Status status;

	int dim, ndata, k;
	double max_double;
	unsigned int seedArray[4];
	unsigned int querySeed;
	if (SEED_SET == 1)
	{
		seedArray[1] = 1;
		seedArray[2] = 2;
		seedArray[3] = 3;
		seedArray[4] = 4;
		querySeed = 722;
	}
	else if (SEED_SET == 2)
	{
		seedArray[1] = 5;
		seedArray[2] = 6;
		seedArray[3] = 7;
		seedArray[4] = 8;
		querySeed = 144;
	}
	else if (SEED_SET == 3)
	{
		seedArray[1] = 9;
		seedArray[2] = 10;
		seedArray[3] = 11;
		seedArray[4] = 12;
		querySeed = 314;
	}
	else
	{
		seedArray[1] = 13;
		seedArray[2] = 14;
		seedArray[3] = 15;
		seedArray[4] = 16;
		querySeed = QSEED;
	}


	printf("\nChecking...\n");

	getCmdArgs(argc, argv, &dim, &ndata, &k,&max_double);
	printf("From world_rank %d, dim = %d, ndata = %d, k = %d, max_double = %lf \n\n", world_rank, dim, ndata, k, max_double);

	int subdomain = ndata / (world_size);
	if (world_rank == world_size - 1)
	{
		subdomain += ndata % (world_size);
	}
	int numSeeds = 4 / world_size;

	printf("From world_rank %d, subdomain = %d numSeeds = %d \n", world_rank, subdomain, numSeeds);

	double * dataArray = (double *)malloc(sizeof(double)*subdomain*dim);
	double * query = (double *)malloc(sizeof(double)*dim);
	//double query[] = { 30.00, 40.00 };

	if (world_rank == 0)
	{
		generateRandomArray(dataArray, subdomain * dim, max_double, numSeeds, &seedArray[world_rank]);
		generateRandomArray(query, dim, max_double,1, &querySeed);
		printf("Query Point in rank = %d =>", world_rank);
		printArrayDoubles(query, 1, dim);
		if (world_size > 1)
		{
			MPI_Bcast(query, 2, MPI_DOUBLE, 0, MCW);
			MPI_Barrier(MCW);
		}
	}
	else
	{
		generateRandomArray(dataArray, subdomain * dim, max_double, numSeeds, &seedArray[world_rank * numSeeds]);
		MPI_Bcast(query, 2, MPI_DOUBLE, 0, MCW);
		MPI_Barrier(MCW);
		printf("Query Point in rank = %d =>",world_rank);
		printArrayDoubles(query, 1, dim);
	}



	//At this point, every process has a data array of the correct size and the query point and all of the arguments.
	//Now begin building the kmeans structure.
/******************************************************************************************************************/

struct kmeans * KM = NULL;

kmeans(&KM,dim,ndata,subdomain,dataArray,k);



	//At this point, every process has a local kmeans struct.
	//The search can now be run.
/*******************************************************************************************************************/

struct stackBase * result = initStack(dim);
int pointsSearched = search(KM,query,result);
printf("%d points searched.\n",pointsSearched);












/**********************************************************************************************************************************/




	//Use brute force search to find the nearest point.

	double * resultA = (double *)malloc(sizeof(double)*dim);
	int i,minLoc;
	for (i = 0; i < dim; i++) { resultA[i] = query[i]; }
	double * LocalBresult = (double *)malloc(sizeof(double)*(dim+1));
	double * allDistPoints = (double *)malloc(sizeof(double)*(dim+1)*world_size);
	double absMinDist;


	LocalBresult[0] = bruteForceSearch(dataArray, query, dim, subdomain, resultA, LocalBresult);

	MPI_Gather(LocalBresult, dim + 1, MPI_DOUBLE, allDistPoints, dim + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (world_rank == 0)
	{
		minLoc = findMinimum(allDistPoints, (dim+1)*world_size, &absMinDist, dim+1);
		printf("minLoc = %d \n",minLoc);
		printArrayDoubles(allDistPoints,world_size,dim+1);
		printf("absMinDist = %lf in rank %d \n",absMinDist,world_rank);
		printArrayDoubles(&allDistPoints[minLoc+1], dim, 1);
		if (checkResult((result->firstNode)->pointArray,&allDistPoints[minLoc + 1],dim))
		{
			printf("THE RESULT IS CORRECT. \n");
		}
		else
		{
			printf("THE RESULT IS NOT CORRECT. \n");
		}
	}





if (dim == 2)
{
	writeResults(dim,ndata, dataArray, KM->cluster_assign);
}


	clearStack(result);
	free(result);
	free(LocalBresult);
	free(allDistPoints);
	free(dataArray);
	free(query);
	MPI_Finalize();
}
