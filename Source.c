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
	generateRandomArray(query, dim, max_double,1, &multSeed,querySeed);
	if (world_rank == 0)
	{
		printf("\n");
		printf("Values dim = %d, ndata = %d, k = %d, max_double = %lf, subdomain = %d numSeeds = %d \n", dim, ndata, k, max_double, subdomain, numSeeds);
		printf("Query Point =>");
		printArrayDoubles(query, 1, dim);
	}
	//if (INIT_K) { printf("Data array: \n"); printDataArray(dataArray,dim,subdomain); }

	//At this point, every process has a data array of the correct size and the query point and all of the arguments.
	//Now begin building the kmeans structure.
/******************************************************************************************************************/

struct kmeans * KM = NULL;

kmeans(&KM,dim,subdomain,dataArray,k,world_rank,world_size);

//if (INIT_K) { displayKM(KM); }

	//At this point, every process has a local kmeans struct.
	//The search can now be run.
/*******************************************************************************************************************/

MPI_Barrier(MCW);
struct stackBase * result = initStack(dim);
int pointsSearched = 0,globPointsSearched = 0;
pointsSearched = search(KM,query,result);
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
	printf("->>>%d points searched in total.\n",globPointsSearched);
}
if (INTERESTING_CASE) printf("->>>>>>>>>>>>>>>%d points searched in total on world_rank %d.\n",pointsSearched,world_rank);

/**********************************************************************************************************************************/




	//Use brute force search to find the nearest point.

	int minLoc,isOneResult;
	double * LocalBresult = (double *)malloc(sizeof(double)*(dim+1));
	double * allDistPoints = (double *)malloc(sizeof(double)*(dim+1)*world_size);
	double absMinDist;

	LocalBresult[0] = bruteForceSearch(dataArray, query, dim, subdomain, LocalBresult);
	MPI_Gather(LocalBresult, dim + 1, MPI_DOUBLE, allDistPoints, dim + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (world_rank == 0)
	{
		struct stackNode * iterator = result->firstNode;
		minLoc = findMinimum(allDistPoints, (dim+1)*world_size, &absMinDist, dim+1);
		printf("Print result stack: \n");
		printStack(result);
		printf("BRUTE FORCE RESULT: \n");
		printf("distance to nearest point = %lf \n",absMinDist);
		printArrayDoubles(&allDistPoints[minLoc+1], dim, 1);
		while (iterator != NULL)
		{
			isOneResult = checkResult(iterator->pointArray,&allDistPoints[minLoc + 1],dim);
			if (isOneResult)
			{
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
	}

if (dim == 2)
{
	write_results(dim,subdomain,KM->data,KM->cluster_assign,k,KM->cluster_centroid,world_rank);
}


	//clearStack(result);
	//free(result);
	free(LocalBresult);
	free(allDistPoints);
	free(dataArray);
	free(query);
	MPI_Finalize();
}
