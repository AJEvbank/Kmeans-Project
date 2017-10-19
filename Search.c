#include "mainHeader.h"

int search(struct kmeans * KM, double * query, struct stackBase * result)
{
  int loop_control = 1,pointCount = 0;
  /* Calculate distances to clusters. */
  double * ClusterDistancesLoc = allocateAndInitializeZeroDouble(KM->k);
  double * ClusterDistancesGlob = allocateAndInitializeZeroDouble(KM->k);
  int nearestCluster = GetClusterDistances(KM,query,ClusterDistancesLoc,ClusterDistancesGlob),candidateCluster;
  int counter = 0;

  /* Iterate: do-while */
  /* Find nearest point in nearest cluster. Distance = straight line - radius */
  do {
      MPI_Barrier(MCW);
      pointCount += GetNearestPoint(KM,result,query,nearestCluster);
      MPI_Barrier(MCW);
      ClusterDistancesGlob[nearestCluster] = INFINITY;
      MPI_Barrier(MCW);
      candidateCluster = GetNearestCluster(ClusterDistancesGlob,KM->k,(result->firstNode)->distance,nearestCluster);
      MPI_Barrier(MCW);
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
      MPI_Barrier(MCW);
      counter++;

  }while(loop_control);
  free(ClusterDistancesLoc);
  free(ClusterDistancesGlob);
  return pointCount;
}

int GetClusterDistances(struct kmeans * KM, double * query, double * ClusterDistancesLoc, double * ClusterDistancesGlob)
{
  int i,nearestCluster,subd = (KM->k)/(KM->world_size),start = (KM->world_rank) * subd;
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
   }

   MPI_Barrier(MCW);
   MPI_Allreduce(
                 ClusterDistancesLoc,
                 ClusterDistancesGlob,
                 (KM->k),
                 MPI_DOUBLE,
                 MPI_MAX,
                 MCW
   );
   MPI_Barrier(MCW);

   for (i = 0; i < KM->k; i++)
   {
     distance = ClusterDistancesGlob[i];
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
  double * topPointLoc = allocateAndInitializeZeroDouble(KM->dim + 2);
  double * pointsGlob = allocateAndInitializeZeroDouble((KM->dim + 2) * KM->world_size);
	int i,j,k, dataPoint,start = (KM->cluster_start)[nearestCluster],size = (KM->cluster_size)[nearestCluster];
	double distanceCalculating;

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
  if (result->stackDepth != 0)
  {
    topPointLoc[0] = (result->firstNode)->distance;
    topPointLoc[1] = (double)(result->firstNode)->cluster;
    for(i = 2,j = 0; i < (KM->dim + 2); i++,j++)
    {
      topPointLoc[i] = ((result->firstNode)->pointArray)[j];
    }
  }
  else
  {
    topPointLoc[0] = INFINITY;
    topPointLoc[1] = -1;
    for(i = 2; i < (KM->dim + 2); i++)
    {
      topPointLoc[i] = -1;
    }
  }

  MPI_Barrier(MCW);
  MPI_Allgather(
                topPointLoc,
                (KM->dim + 2),
                MPI_DOUBLE,
                pointsGlob,
                (KM->dim + 2),
                MPI_DOUBLE,
                MCW
  );
  MPI_Barrier(MCW);

  SetStackWithGlobal(KM,result,pointsGlob);

	free(thisPoint);
  free(topPointLoc);
  free(pointsGlob);
	return (KM->cluster_size)[nearestCluster];
}

int GetNearestCluster(double * ClusterDistancesGlob, int size, double minDist, int nearestCluster)
{
  int i;
  for (i = 0; i < size; i++)
  {
    if (ClusterDistancesGlob[i] <= minDist)
    {
      return i;
    }
  }
  return nearestCluster;
}

void SetStackWithGlobal(struct kmeans * KM, struct stackBase * result, double * pointsGlob)
{
  int i,first_index,bestCluster;
  double minDist;
  if (result->stackDepth == 0)
  {
    minDist = INFINITY;
  }
  else
  {
    minDist = (result->firstNode)->distance;
  }
  struct stackNode * iterator = NULL;
  for (i = 0; i < KM->world_size; i++)
  {
    first_index = i * (KM->dim + 2);
    if (pointsGlob[first_index] < minDist)
    {
      minDist = pointsGlob[first_index];
      bestCluster = pointsGlob[first_index + 1];
      clearStack(result);
      pushNode(&pointsGlob[first_index + 2],minDist,bestCluster,result);
    }
    else if (pointsGlob[first_index] == minDist)
    {
      minDist = pointsGlob[first_index];
      bestCluster = pointsGlob[first_index + 1];
      iterator = result->firstNode;
      while (iterator != NULL)
      {
        if(checkResult(iterator->pointArray,&pointsGlob[first_index + 2],KM->dim))
        {
          break;
        }
        iterator = iterator->nextNode;
      }
      if (iterator == NULL)
      {
          pushNode(&pointsGlob[first_index + 2],minDist,bestCluster,result);
      }
    }
  }
  return;
}
