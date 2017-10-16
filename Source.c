#include "mainHeader.h"

int main(int argc, char** argv) {

	// // Initialize the MPI environment
	// MPI_Init(&argc, &argv);
	// // Find out rank, size
	// int world_rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// int world_size;
	// MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	// //MPI_Status status;

	int dim, ndata, k;
	double max_double;
	unsigned int seedArray[] = { 2, 4, 8, 16 };
	unsigned int querySeed;
	int seedMult,numSeeds = 4;


	getCmdArgs(argc, argv, &dim, &ndata, &k,&max_double,&seedMult);

	printf("Values: dim = %d, ndata = %d, k = %d, max_double = %lf, ndata = %d numSeeds = %d, seedMult = %d  \n\n", dim, ndata, k, max_double, ndata, numSeeds,seedMult);

	double * dataArray = (double *)malloc(sizeof(double)*ndata*dim);
	double * query = (double *)malloc(sizeof(double)*dim);
	generateRandomArray(dataArray, ndata * dim, max_double, numSeeds, &seedArray[0]);
	querySeed = seedMult;
	generateRandomArray(query, dim, max_double,1, &querySeed);
	printf("Query Point: \n");
	printArrayDoubles(query, 1, dim);


	//At this point, every process has a data array of the correct size and the query point and all of the arguments.
	//Now begin building the kmeans structure.
/******************************************************************************************************************/

struct kmeans * KM = NULL;

kmeans(&KM,dim,ndata,ndata,dataArray,k);


	//At this point, every process has a local kmeans struct.
	//The search can now be run.
/*******************************************************************************************************************/

struct stackBase * result = initStack(dim);
int pointsSearched = search(KM,query,result);
printf("->>>%d points searched.\n",pointsSearched);


/**********************************************************************************************************************************/




	//Use brute force search to find the nearest point.


	int isOneResult;
	double * LocalBresult = (double *)malloc(sizeof(double)*(dim+1));


	LocalBresult[0] = bruteForceSearch(dataArray, query, dim, ndata, LocalBresult);

	struct stackNode * iterator = result->firstNode;
	printf("\n");
	printf("Print result stack from kmeans search function: \n");
	printStack(result);
	printf("BRUTE FORCE RESULT: \n");
	printArrayDoubles(&LocalBresult[1], dim, 1);
	printf("minimum distance = %lf \n",LocalBresult[0]);
	while (iterator != NULL)
	{
		isOneResult = checkResult(iterator->pointArray,&LocalBresult[1],dim);
		if (isOneResult)
		{
			break;
		}
		else
		{
			iterator = iterator->nextNode;
		}
	}
	printf("\n");
	if (isOneResult)
	{
		printf("THE RESULT IS CORRECT. \n");
	}
	else
	{
		printf("THE RESULT IS NOT CORRECT. \n");
	}
	printf("\n");

if (dim == 2)
{
	writeResults(dim,ndata, dataArray, KM->cluster_assign);
}


	clearStack(result);
	free(result);
	free(LocalBresult);
	free(dataArray);
	free(query);
	//MPI_Finalize();
}
