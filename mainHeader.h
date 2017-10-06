#ifndef MAINHEADER

#define MAINHEADER


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "math.h"


/* MACROS */

#define DIM 2
#define NDATA 40
#define K 4
#define MAX_DOUBLE 50.00
#define MCW MPI_COMM_WORLD
#define QSEED 30
#define SEED_SET 1
#define FIRST_CENTROID_SEED 13
#define THRESHOLD 25




#define DEBUG 0
#define DEBUG_RANDOM 0
#define DISPLAY_KM_INIT 1
#define DEBUG_SELECTK 0
#define DEBUG_ASSIGN 1
#define DEBUG_THRESHOLD 1

#define WRITE_RESULTS 1

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
	int * cluster_size;
	int * cluster_start;
	double * cluster_radius;
	double ** cluster_centroid;
	int * cluster_assign;
	int * cluster_group;
};

/* Prototypes */

void getCmdArgs(int argc, char ** argv, int * dim, int * ndata, int * k, double * max_double);

int isNumber(const char * str);

void printArray(int * nums, int count);

void printArrayDoubles(double * nums, int ndata, int dim);

int generateRandomArray(double * dataArray, int subdomain, double max_double, int seeds, unsigned int * seedArray);

double bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * result, double * Bresult);

int findMinimum(double * Array, int size, double * minimum, int stride);

/* initialize the kmeans structure */

void initializeKM(struct kmeans ** KM, int dim, int ndata, double * dataArray, int k);

int * allocateAndInitializeZeroInt(int size_of_target);

double ** allocateAndInitializeZeroDoubleMulti(int k, int dimension);

double * allocateAndInitializeZeroDouble(int size_of_target);

/* debugging displays */

void displayKM(struct kmeans * KM);

void printArraysInt(int * ArrayInt, int size, const char * text);

void printArraysDouble(double ** ArrayDouble, int size, int dim, const char * text);

void printArrayKMD(double * array, int size);

void printDataArray(double * dataArray, int dim, int ndata);

void writeResults(int dim, int ndata, double* data, int* cluster_assign);

/* GetKCentroids.c */

void GetKCentroids(struct kmeans * KM);

int GetNextCluster(struct kmeans * KM, int numClusters);

double GetDistance2PointsDC(struct kmeans *KM, int first_index, int centroid);

/* ClusterizeKM.c */

void ClusterizeKM(struct kmeans * KM, int threshold);

void AssignDPs(struct kmeans * KM);

void RecalculateCentroids(struct kmeans * KM);


#endif
