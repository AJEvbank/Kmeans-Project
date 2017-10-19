#ifndef MAINHEADER

#define MAINHEADER


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "math.h"


/* MACROS */

#define DIM 128
#define NDATA 10000
#define K 128
#define MAX_DOUBLE 1000.00
#define MCW MPI_COMM_WORLD
#define SEEDMULT 1
#define QSEED 666

#define INTERESTING_CASE 0
#define INIT_K 1

/* Data Structures */

enum isNumStates {
	INITIAL,
	PLUSORMINUS,
	ZERO,
	NUMBER,
	DECIMAL,
	ERROR
};

struct kmeans {
	int dim;
	int ndata;
	double * data;
	int k;
	int world_rank;
	int world_size;
	int * cluster_size;
	int * cluster_start;
	double * cluster_radius;
	double ** cluster_centroid;
	int * cluster_assign;
	//int * cluster_group;
};

struct stackBase {
	int stackDepth;
	int arraySize;
	struct stackNode * firstNode;
};

struct stackNode {
	double * pointArray;
	double distance;
	int cluster;
	struct stackNode * nextNode;
};

/* Prototypes */

void getCmdArgs(int argc, char ** argv, int * dim, int * ndata, int * k, double * max_double, unsigned int * multSeed);

int isNumber(const char * str);

void printArray(int * nums, int count);

void printArrayDoubles(double * nums, int ndata, int dim);

void generateRandomArray(double * dataArray, int domain, double max_double, int seeds, unsigned int * seedArray, unsigned int seedMult);

double bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * Bresult);

int findMinimum(double * Array, int size, double * minimum, int stride);

int checkResult(double * searchResult, double * bruteResult, int dim);

/* initialize the kmeans structure */

void initializeKM(struct kmeans ** KM, int dim, int ndata, double * dataArray, int k, int world_rank, int world_size);

int * allocateAndInitializeZeroInt(int size_of_target);

double ** allocateAndInitializeZeroDoubleMulti(int k, int dimension);

double * allocateAndInitializeZeroDouble(int size_of_target);

void kmeans(struct kmeans ** KM, int dim, int ndata, double * dataArray, int k, int world_rank, int world_size);

/* debugging displays - DEBUG.c */

#define SHOW_DP_NUMBER 1

void displayKM(struct kmeans * KM);

void printArraysInt(int * ArrayInt, int size, const char * text);

void printArraysDouble(double ** ArrayDouble, int size, int dim, const char * text);

void printArrayKMD(double * array, int size);

void printDataArray(double * dataArray, int dim, int ndata);

void write_results(int dim, int ndata, double *data, int *cluster_assign, int k, double **cluster_centroids, int world_rank);

void printArrayDouble(double * ArrayDouble, int size, const char * text);

void printStack(struct stackBase *stack);

void displayAverageDistance(double * clusterDistances, int size);

/* GetKCentroids.c */

void GetKCentroids(struct kmeans * KM);

int GetNextCluster(struct kmeans * KM, int numClusters);

double GetDistance2PointsDC(struct kmeans *KM, int first_index, int centroid);

/* ClusterizeKM.c */

void ClusterizeKM(struct kmeans * KM);

void AssignDPs(struct kmeans * KM);

int RecalculateCentroids(struct kmeans * KM);

void SaveClusters(struct kmeans * KM);


/* QuickSortCode.c */

void quickSort( struct kmeans * KM, int l, int r);

int partition( struct kmeans * KM, int l, int r);

void SwapPoints(struct kmeans * KM,int left,int right);

void SortDataArray(struct kmeans * KM);

void CalculateRadii(struct kmeans * KM);

void EmptyClusters(struct kmeans * KM);

void DeleteEmptyCluster(struct kmeans * KM, int i);

/* Search.c */

int search(struct kmeans * KM, double * query, struct stackBase * result);

int GetClusterDistances(struct kmeans * KM, double * query, double * ClusterDistancesLoc, double * ClusterDistancesGlob);

double GetDistance2PointsQC(struct kmeans * KM, double * query, int cluster);

int GetNearestPoint(struct kmeans * KM, struct stackBase * result, double * query, int nearestCluster);

int GetNearestCluster(double * ClusterDistances, int size, double minDist, int nearestCluster);

void SetStackWithGlobal(struct kmeans * KM, struct stackBase * result, double * pointsGlob);

/* Stack.c */

struct stackBase * initStack(int);

void pushNode(double *, double, int,  struct stackBase *);

void pop(struct stackBase *);

int peekDepth(struct stackBase * );

void clearStack(struct stackBase *);

#endif
