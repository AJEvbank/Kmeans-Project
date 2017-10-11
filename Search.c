#include "mainHeader.h"

int search(struct kmeans * KM, double * query, struct stackBase * result)
{
  int loop_control = 1,pointCount = 0;
  /* Calculate distances to clusters. */
  double * ClusterDistancesLoc = allocateAndInitializeZeroDouble(KM->k);
  int nearestCluster = GetClusterDistances(KM,query,ClusterDistancesLoc),candidateCluster;
  if (DEBUG_CLUSTER_DIST) { printArrayDouble(ClusterDistances,KM->k,"distance to cluster"); }

  /* Iterate: do-while */
  /* Find nearest point in nearest cluster. Distance = straight line - radius */
  do {
      pointCount += GetNearestPoint(KM,result,query,nearestCluster);
      if (QUERY_ANALYSIS) { printStack(result); }
      ClusterDistances[nearestCluster] = INFINITY;
      candidateCluster = GetNearestCluster(ClusterDistances,KM->k,(result->firstNode)->distance,nearestCluster);
      if (DEBUG_CLUSTER_DIST) { printArrayDouble(ClusterDistances,KM->k,"distance to cluster"); }
      if (DEBUG_CLUSTER_DIST) { printf("candidateCluster %d \n",candidateCluster); }
      /* If there is a cluster K nearer than the current nearest point,
         then reiterate searching K to update the nearest point. */
      if (candidateCluster != nearestCluster)
      {
        nearestCluster = candidateCluster;
        if (DEBUG_CLUSTER_DIST) { printf("nearestCluster = %d, loop_control = %d \n",nearestCluster,loop_control); }
      }
      else
      {
        loop_control = 0;
      }
  }while(loop_control);
  free(ClusterDistancesLoc);
  return pointCount;
}

int GetClusterDistances(struct kmeans * KM, double * query, double * ClusterDistancesLoc)
{
  int i,nearestCluster,subd = (KM->k)/world_size,start = world_rank * subd;
  double distance, minDist = INFINITY;
  for (i = start; i < (subd + start); i++)
  {
    distance = GetDistance2PointsQC(KM,query,i);
    distance = distance - (KM->cluster_radius)[i];
    if (distance < 0)
    {
      distance = 0;
    }
    ClusterDistancesLoc[i] = distance;
    if (distance < minDist)
    {
      minDist = distance;
      nearestCluster = i;
    }
  }
  return nearestCluster;
}

double GetDistance2PointsQC(struct kmeans * KM, double * query, int cluster)
{
  int i;
  double distCalc = 0;
  for (i = 0; i < KM->dim; i++)
  {
    distCalc += pow( fabs(query[i] - (KM->cluster_centroid)[cluster][i]), 2);
  }
  distCalc = sqrt(distCalc);
  return distCalc;
}

int GetNearestPoint(struct kmeans * KM, struct stackBase * result, double * query, int nearestCluster)
{
  double * thisPoint = (double *)malloc(sizeof(double)*KM->dim);
	int i,j,k, dataPoint,start = (KM->cluster_start)[nearestCluster],size = (KM->cluster_size)[nearestCluster];
	double distanceCalculating;
	// if (QUERY_ANALYSIS2) { printf("(Tree->cluster_start)[nearestCluster] = %d,\n", start);
	// printf("(Tree->cluster_size)[nearestCluster] = %d,\n", size);
	// }
	for (i = start; i < start + size; i++)
	{
		dataPoint = i * KM->dim;
		distanceCalculating = 0;
		// if (QUERY_ANALYSIS2) { printf("dataPoint = %d \n", dataPoint); }
		for (j = 0; j < KM->dim; j++)
		{
			distanceCalculating += pow( fabs((query[j] - (KM->data)[dataPoint+j])), 2);
			// if (QUERY_ANALYSIS2) { printf("with data[%d] -> %lf \n",(dataPoint+j), distanceCalculating); }
		}
		distanceCalculating = sqrt(distanceCalculating);
		// if (QUERY_ANALYSIS2) { printf("final = %lf \n", distanceCalculating); }
		if (result->stackDepth == 0) {
			for (k = 0 ; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
			// if (QUERY_ANALYSIS) { printArray(thisPoint, Tree->dim); }
			pushNode(thisPoint, distanceCalculating, nearestCluster, result);
			// if (QUERY_ANALYSIS) { printStack(result); }
		}
		else
		{
			if (distanceCalculating < (result->firstNode)->distance)
			{
				for (k = 0; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
				// if (QUERY_ANALYSIS) { printArray(thisPoint, Tree->dim); }
				clearStack(result);
				pushNode(thisPoint, distanceCalculating, nearestCluster, result);
			}
			else if (distanceCalculating == (result->firstNode)->distance)
			{
				for (k = 0; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
				// if (QUERY_ANALYSIS) { printArray(thisPoint, Tree->dim); }
				pushNode(thisPoint, distanceCalculating, nearestCluster, result);
			}
		}
	}
	// if (QUERY_ANALYSIS) { printStack(result); }
	free(thisPoint);
	return (KM->cluster_size)[nearestCluster];
}

int GetNearestCluster(double * ClusterDistances, int size, double minDist, int nearestCluster)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (ClusterDistances[i] <= minDist)
    {
      return i;
    }
  }
  return nearestCluster;
}
