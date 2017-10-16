#include "mainHeader.h"

int search(struct kmeans * KM, double * query, struct stackBase * result)
{
  int loop_control = 1,pointCount = 0;
  /* Calculate distances to clusters. */
  double * ClusterDistances = allocateAndInitializeZeroDouble(KM->k);
  int nearestCluster = GetClusterDistances(KM,query,ClusterDistances),candidateCluster;

  /* Iterate: do-while */
  /* Find nearest point in nearest cluster. Distance = straight line - radius */
  do {
      pointCount += GetNearestPoint(KM,result,query,nearestCluster);
      ClusterDistances[nearestCluster] = INFINITY;
      candidateCluster = GetNearestCluster(ClusterDistances,KM->k,(result->firstNode)->distance,nearestCluster);
      /* If there is a cluster K nearer than the current nearest point,
         then reiterate searching K to update the nearest point. */
      if (candidateCluster != nearestCluster)
      {
        nearestCluster = candidateCluster;
      }
      else
      {
        loop_control = 0;
      }
  }while(loop_control);
  return pointCount;
}

int GetClusterDistances(struct kmeans * KM, double * query, double * ClusterDistances)
{
  int i,nearestCluster;
  double distance, minDist = INFINITY;
  for (i = 0; i < KM->k; i++)
  {
    distance = GetDistance2PointsQC(KM,query,i);
    distance = distance - (KM->cluster_radius)[i];
    if (distance < 0)
    {
      distance = 0;
    }
    ClusterDistances[i] = distance;
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
	// }
	for (i = start; i < start + size; i++)
	{
		dataPoint = i * KM->dim;
		distanceCalculating = 0;
		for (j = 0; j < KM->dim; j++)
		{
			distanceCalculating += pow( fabs((query[j] - (KM->data)[dataPoint+j])), 2);
		}
		distanceCalculating = sqrt(distanceCalculating);
		if (result->stackDepth == 0) {
			for (k = 0 ; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
			pushNode(thisPoint, distanceCalculating, nearestCluster, result);
		}
		else
		{
			if (distanceCalculating < (result->firstNode)->distance)
			{
				for (k = 0; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
				clearStack(result);
				pushNode(thisPoint, distanceCalculating, nearestCluster, result);
			}
			else if (distanceCalculating == (result->firstNode)->distance)
			{
				for (k = 0; k < KM->dim; k++) { thisPoint[k] = (KM->data)[dataPoint + k]; }
				pushNode(thisPoint, distanceCalculating, nearestCluster, result);
			}
		}
	}
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
