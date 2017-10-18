#ifndef COMMANDLINEARGS

#define COMMANDLINEARGS


#include "mainHeader.h"

void getCmdArgs(int argc, char ** argv, int * dim, int * ndata, int * k, double * max_double,int * seedMult)
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
	if (argc >= 6)
	{
		if (isNumber(argv[5]))
		{
			*seedMult = (int)atof(argv[5]);
		}
		if (*seedMult <= 0)
		{
			*seedMult = SEEDMULT;
		}
	}
	else
	{
		*seedMult = SEEDMULT;
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
		for (j = 0; j < domain / seeds; j++)
		{
			dataArray[first_index + j] = ((double)rand() / (double)RAND_MAX) * max_double;
		}
	}
	return 0;
}

double bruteForceSearch(double * dataArray, double * query, int dim, int ndata, double * Bresult)
{
	int first_index, i, j, nearestPoint;
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
	}
	return minDist;
}

int checkResult(double * searchResult, double * bruteResult, int dim)
{
	int i;
	for (i = 0; i < dim; i++)
	{
		if(searchResult[i] != bruteResult[i])
		{
			return 0;
		}
	}
	return 1;
}

void setSeedArray(unsigned int * seedArray, int seedMult)
{
	int i;
	for (i = 0; i < 4; i++)
	{
		seedArray[i] = seedArray[i] * seedMult;
	}
	return;
}



#endif
