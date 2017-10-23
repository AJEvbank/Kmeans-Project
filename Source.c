#include "mainHeader.h"

int main(int argc, char** argv) {

	// Set up clock check
	clock_t begin_build_kmeans;
	clock_t end_build_kmeans;
	clock_t begin_search_kmeans;
	clock_t end_search_kmeans;
	clock_t begin_brute_force_search;
	clock_t end_brute_force_search;
	clock_t global_minimum_start;
	clock_t global_maximum_end;

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
	unsigned int seedArray[] = { 2, 4, 8, 16 };
	unsigned int multSeed,querySeed = QSEED;

	getCmdArgs(argc, argv, &dim, &ndata, &k,&max_double,&multSeed);

	int subdomain = ndata / (world_size);
	if (world_rank == world_size - 1)
	{
		subdomain += ndata % (world_size);
	}
	int numSeeds = 4 / world_size;

	double * dataArray = (double *)malloc(sizeof(double)*subdomain*dim);
	double * query = (double *)malloc(sizeof(double)*dim);

	generateRandomArray(dataArray, subdomain * dim, max_double, numSeeds, &seedArray[world_rank * numSeeds], multSeed);
	generateRandomArray(query, dim, max_double,1, &querySeed,multSeed);
	if (world_rank == 0)
	{
		printf("\n");
		printf("Values dim = %d, ndata = %d, k = %d, max_double = %lf, subdomain = %d, numSeeds = %d, cycles per second = %ld \n\n", dim, ndata, k, max_double, subdomain, numSeeds,CLOCKS_PER_SEC);
		printf("Query Point =>");
		printArrayDoubles(query, 1, dim);
	}

	//At this point, every process has a data array of the correct size and the query point and all of the arguments.
	//Now begin building the kmeans structure.
/******************************************************************************************************************/

begin_build_kmeans = clock();

struct kmeans * KM = NULL;

kmeans(&KM,dim,subdomain,dataArray,k,world_rank,world_size);

end_build_kmeans = clock();
MPI_Barrier(MCW);

MPI_Reduce(
						&begin_build_kmeans,
						&global_minimum_start,
						1,
						MPI_LONG,
						MPI_MIN,
						0,
						MCW
);
MPI_Barrier(MCW);
MPI_Reduce(
						&end_build_kmeans,
						&global_maximum_end,
						1,
						MPI_LONG,
						MPI_MAX,
						0,
						MCW
);
MPI_Barrier(MCW);

if (world_rank == 0)
{
	printf("\n");
	printf("begin build = %ld cycles \n",global_minimum_start);
	printf("end build = %ld cycles \n",global_maximum_end);
	printf("Cycles to build kmeans structure = %ld \n",(global_maximum_end - global_minimum_start));
	printf("Time to build kmeans structure = %lf seconds \n",((double)(global_maximum_end - global_minimum_start))/(double)CLOCKS_PER_SEC);
}

if (KM_TEST1) displaySelectedFromKM(KM,1,0,1,1,1,1,0);
	//At this point, every process has a local kmeans struct.
	//The search can now be run.
/*******************************************************************************************************************/

MPI_Barrier(MCW);
begin_search_kmeans = clock();
global_minimum_start = begin_search_kmeans;

struct stackBase * result = initStack(dim);
int pointsSearched = 0,globPointsSearched = 0;
pointsSearched = search(KM,query,result);

do {
	end_search_kmeans = clock();
}while(end_search_kmeans == begin_search_kmeans);
global_maximum_end = end_search_kmeans;
MPI_Barrier(MCW);

MPI_Reduce(
						&begin_search_kmeans,
						&global_minimum_start,
						1,
						MPI_LONG,
						MPI_MIN,
						0,
						MCW
);
MPI_Barrier(MCW);
MPI_Reduce(
						&end_search_kmeans,
						&global_maximum_end,
						1,
						MPI_LONG,
						MPI_MAX,
						0,
						MCW
);
MPI_Barrier(MCW);

MPI_Reduce(
						&pointsSearched,
						&globPointsSearched,
						1,
						MPI_INT,
						MPI_SUM,
						0,
						MCW
);
MPI_Barrier(MCW);
if (world_rank == 0)
{
	printf("\n");
	printf("begin search = %ld cycles \n",global_minimum_start);
	printf("end search = %ld cycles \n",global_maximum_end);
	printf("Cycles to search kmeans structure = %ld \n",(global_maximum_end - global_minimum_start));
	printf("Time to search kmeans structure = %lf seconds \n",((double)(global_maximum_end - global_minimum_start))/(double)CLOCKS_PER_SEC);
	printf("->>>%d points searched in total.\n",globPointsSearched);
}

/**********************************************************************************************************************************/

	//Use brute force search to find the nearest point.

	MPI_Barrier(MCW);
	begin_brute_force_search = clock();
	global_minimum_start = begin_brute_force_search;

	int minLoc,isOneResult;
	double * LocalBresult = (double *)malloc(sizeof(double)*(dim+1));
	double * allDistPoints = (double *)malloc(sizeof(double)*(dim+1)*world_size);
	double absMinDist;

	LocalBresult[0] = bruteForceSearch(dataArray, query, dim, subdomain, LocalBresult);
	MPI_Gather(
							LocalBresult,
							dim + 1,
							MPI_DOUBLE,
							allDistPoints,
							dim + 1,
							MPI_DOUBLE,
							0,
							MPI_COMM_WORLD
	);

	do {
		end_brute_force_search = clock();
	}while(end_brute_force_search == begin_brute_force_search);
	global_maximum_end = end_brute_force_search;
	MPI_Barrier(MCW);

	MPI_Reduce(
							&begin_brute_force_search,
							&global_minimum_start,
							1,
							MPI_LONG,
							MPI_MIN,
							0,
							MCW
	);
	MPI_Barrier(MCW);
	MPI_Reduce(
							&end_brute_force_search,
							&global_maximum_end,
							1,
							MPI_LONG,
							MPI_MAX,
							0,
							MCW
	);
	MPI_Barrier(MCW);

	if (world_rank == 0)
	{
		printf("\n");
		printf("begin brute force search = %ld cycles \n",global_minimum_start);
		printf("end brute force search = %ld cycles \n",global_maximum_end);
		printf("Cycles to search kmeans structure = %ld \n",(global_maximum_end - global_minimum_start));
		printf("Time for brute force search = %lf seconds \n\n",((double)(global_maximum_end - global_minimum_start))/(double)CLOCKS_PER_SEC);
	}


	if (world_rank == 0)
	{
		int i;
		struct stackNode * iterator = result->firstNode;
		minLoc = findMinimum(allDistPoints, (dim+1)*world_size, &absMinDist, dim+1);

		printf("Result stack date: \n");
		printf("stackDepth = %d, firstNode = %p \n", result->stackDepth,(void *)result->firstNode);
		printf("distance to query point = %lf \n", result->firstNode->distance);

		printf("\n");

		printf("BRUTE FORCE RESULT: \n");
		printf("distance to nearest point = %lf \n",absMinDist);

		while (iterator != NULL)
		{
			isOneResult = checkResult(iterator->pointArray,&allDistPoints[minLoc + 1],dim);
			if (isOneResult)
			{
				printf("\n");
				printf("Search Result: \t\t Brute Force Result: \n");
				for (i = 0; i < dim; i++)
				{
					printf("%lf \t\t %lf \n",(iterator->pointArray)[i],allDistPoints[minLoc + 1 + i]);
				}
				printf("In cluster %d \n",iterator->cluster);
				printf("\n");
				break;
			}
			else
			{
				iterator = iterator->nextNode;
			}
		}
		if (isOneResult)
		{
			printf("THE RESULT IS CORRECT. \n");
		}
		else
		{
			printf("THE RESULT IS NOT CORRECT. \n");
		}
		printf("\n");
		if (result->stackDepth > 1)
		{
			printf("Kmeans search stack depth greater than 1. \n");
			printf("Stack: \n");
			printStack(result);
		}
	}

if (dim == 2)
{
	write_results(dim,subdomain,KM->data,KM->cluster_assign,k,KM->cluster_centroid,world_rank);
}


	clearStack(result);
	free(result);
	free(LocalBresult);
	free(allDistPoints);
	free(dataArray);
	free(query);
	MPI_Finalize();
}
