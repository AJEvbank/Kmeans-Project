#include "mainHeader.h"

void GetKCentroids(struct kmeans * KM)
{
  /* 1. Randomly select a point in the data array for the first centroid. */
  int first_centroid = 0;
  int i;
  for (i = 0; i < KM->dim; i++)
  {
    (KM->cluster_centroid)[0][i] = (KM->data)[(first_centroid * (KM->dim)) + i];
  }


  /* 2. Select each subsequent centroid based on the max-of-the-min criteria given in class. */




  return;
}
