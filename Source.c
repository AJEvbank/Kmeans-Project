#include "mainHeader.h"

int main(int argc, const char** argv) {

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Find out rank, size
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Status status;

	int dim, ndata, k;
	double max_double;
	unsigned int seedArray[4] = { 1, 2, 3, 4 };
	unsigned int querySeed = QSEED;

	printf("\nChecking...\n");

	if (world_rank == 0)
	{
		getCmdArgs(argc, argv, &dim, &ndata, &k,&max_double);
		if (world_size > 1)
		{
			MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&ndata, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&max_double, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		printf("From world_rank %d, dim = %d, ndata = %d, k = %d max_double = %lf \n\n", world_rank, dim, ndata, k, max_double);
	}
	else
	{
		MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&ndata, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&max_double, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		printf("From world_rank %d, dim = %d, ndata = %d, k = %d max_double = %lf \n\n", world_rank, dim, ndata, k, max_double);
	}

	int subdomain = ndata / (world_size);
	if (world_rank == world_size - 1)
	{
		subdomain += ndata % (world_size);
	}
	int numSeeds = 4 / world_size;

	printf("From world_rank %d, subdomain = %d numSeeds = %d \n", world_rank, subdomain, numSeeds);

	double * dataArray = (double *)malloc(sizeof(double)*subdomain*dim);
	double * query = (double *)malloc(sizeof(double)*dim);
	
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
	
	if (DEBUG_RANDOM)
	{
		int i, j, first_index;
		printf("WORLD_RANK = %d \n", world_rank);
		printf("[ %d \n", world_rank);
		for (i = 0; i < subdomain; i++)
		{
			first_index = dim * i;
			for (j = 0; j < dim; j++)
			{
				printf("%lf, ", dataArray[first_index + j]);
			}
			printf("\t%d\n", world_rank);
		}
		printf("] %d \n", world_rank);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//At this point, every process has a data array of the correct size and the query point and all of the arguments.
	//Use brute force search to find the nearest point.

	double * result = (double *)malloc(sizeof(double)*dim);
	int i,minLoc;
	for (i = 0; i < dim; i++) { result[i] = query[i]; }
	double * Bresult = (double *)malloc(sizeof(double)*dim);
	double * minDistances = (double *)malloc(sizeof(double)*world_size);

	double minDist = bruteForceSearch(dataArray, query, dim, subdomain, result, Bresult);
	double absMinDist = 0;

		
	MPI_Allgather(&minDist, 1, MPI_DOUBLE, minDistances, 1, MPI_DOUBLE, MCW);

	if (world_rank == 0)
	{
		minLoc = findMinimum(minDistances, world_size, &absMinDist);
	}
	MPI_Bcast(&absMinDist, 1, MPI_DOUBLE, 0, MCW);
	MPI_Barrier(MCW);

	

	printf("minDistances in rank %d: \n",world_rank);
	printArrayDoubles(minDistances, world_size, 1);
	printf("absMinDist = %lf in rank %d \n",absMinDist,world_rank);
	
	

	


	
	
	

	MPI_Finalize();
}