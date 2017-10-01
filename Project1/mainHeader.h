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
#define DEBUG 0
#define DEBUG_RANDOM 1

/* Data Structures */

enum isNumStates {
	INITIAL,
	PLUSORMINUS,
	ZERO,
	NUMBER,
	DECIMAL,
	ERROR
};

/* Prototypes */

void getCmdArgs(int argc, const char ** argv, int * dim, int * ndata, int * k, double * max_double);

int isNumber(const char * str);

void printArray(int * nums, int count);

void printArrayDoubles(double * nums, int ndata, int dim);

int generateRandomArray(double * dataArray, int subdomain, double max_double, int seeds, int * seedArray);

double bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * result, double * Bresult);

int findMinimum(double * Array, int size, double * minimum);

#endif