#ifndef COMMANDLINEARGS

#define COMMANDLINEARGS


#include "mainHeader.h"

void getCmdArgs(int argc, char ** argv, int * dim, int * ndata, int * k, double * max_double)
{
	if (argc >= 2)
	{
		if (isNumber(argv[1]))
		{
			*dim = (int)atof(argv[1]);
		}
		if (*dim <= 0)
		{
			*dim = DIM;
		}
	}
	else
	{
		*dim = DIM;
	}
	if (argc >= 3)
	{
		if (isNumber(argv[2]))
		{
			*ndata = (int)atof(argv[2]);
		}
		if (*ndata <= 0)
		{
			*ndata = NDATA;
		}
	}
	else
	{
		*ndata = NDATA;
	}
	if (argc >= 4)
	{
		if (isNumber(argv[3]))
		{
			*max_double = atof(argv[3]);
		}
		if (*max_double <= 0)
		{
			*max_double = MAX_DOUBLE;
		}
	}
	else
	{
		*max_double = MAX_DOUBLE;
	}
	if (argc >= 5)
	{
		if (isNumber(argv[4]))
		{
			*k = (int)atof(argv[4]);
		}
		if (*k <= 0)
		{
			*k = K;
		}
	}
	else
	{
		*k = K;
	}

	return;
}

int isNumber(const char * str)
{
	int i;
	int len = strlen(str);
	char ch;
	enum isNumStates state = INITIAL;
	if (len > 0)
	{
		for (i = 0; i < len; i++)
		{
			ch = str[i];
			switch (state)
			{
			case INITIAL:
			{
				if (ch == '-' || ch == '+')
				{
					state = PLUSORMINUS;
				}
				else if (ch == '0')
				{
					state = ZERO;
				}
				else if (isdigit(ch))
				{
					state = NUMBER;
				}
				else
				{
					state = ERROR;
				}
				break;
			}

			case PLUSORMINUS:
			{
				if (ch == '0')
				{
					state = ZERO;
				}
				else if (isdigit(ch))
				{
					state = NUMBER;
				}
				else
				{
					state = ERROR;
				}
				break;
			}

			case ZERO:
			{
				if (ch == '.')
				{
					state = DECIMAL;
				}
				else
				{
					state = ERROR;
				}
				break;
			}

			case NUMBER:
			{
				if (isdigit(ch))
				{
					state = NUMBER;
				}
				else if (ch == '.')
				{
					state = DECIMAL;
				}
				else
				{
					state = ERROR;
				}
				break;
			}

			case DECIMAL:
			{
				if (isdigit(ch))
				{
					state = DECIMAL;
				}
				else
				{
					state = ERROR;
				}
				break;
			}

			case ERROR:
			{
				return 0;
				break;
			}

			default:
			{
				printf("default in isNumber: %d \n", (int)state);
				return 0;
				break;
			}

			}
		}
		return 1;
	}
	else
	{
		return 0;
	}
}

int generateRandomArray(double * dataArray, int domain, double max_double, int seeds, unsigned int * seedArray)
{
	int i,j,first_index;
	for (i = 0; i < seeds; i++)
	{
		srand(seedArray[i]);
		first_index = i * domain/seeds;
		if (DEBUG) { printf("seed = %d \n", seedArray[i]); }
		for (j = 0; j < domain / seeds; j++)
		{
			dataArray[first_index + j] = ((double)rand() / (double)RAND_MAX) * max_double;
		}
	}
	return 0;
}

double bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * result, double * Bresult)
{
	int first_index, i, j, nearestPoint;
	//int isResult = 1;
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
	first_index = nearestPoint * dim;
	for (i = 0; i < dim; i++)
	{
		Bresult[i+1] = dataArray[(first_index) + i];
		// if (dataArray[(first_index) + i] != result[i])
		// {
		// 	isResult = 0;
		// }
	}
	return minDist;
}

int findMinimum(double * Array, int size, double * minimum, int stride)
{
	int minIndex = 0, i;
	*minimum = INFINITY;
	for (i = 0; i < size; i+=stride)
	{
		printf("Array[%d] = %lf \n",i,Array[i]);
		if (Array[i] < *minimum)
		{
			*minimum = Array[i];
			minIndex = i;
		}
	}
	return minIndex;
}







#endif
