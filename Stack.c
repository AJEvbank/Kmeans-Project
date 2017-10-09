#include "mainHeader.h"

struct stackBase * initStack(int size)
{
	struct stackBase * rtrn = (struct stackBase *)malloc(sizeof(struct stackBase));
	rtrn->firstNode = NULL;
	rtrn->stackDepth = 0;
	rtrn->arraySize = size;
	return rtrn;
}

void pushNode(double * rPoint, double dist, int cluster, struct stackBase * thisSB)
{
	int i;
	struct stackNode * temp = (struct stackNode *)malloc(sizeof(struct stackNode));
	double * pointArray = (double *)malloc(sizeof(double)*(thisSB->arraySize));
	for (i = 0; i < thisSB->arraySize; i++) { pointArray[i] = rPoint[i]; }
	temp->nextNode = thisSB->firstNode;
	temp->pointArray = pointArray;
	temp->distance = dist;
	temp->cluster = cluster;
	thisSB->firstNode = temp;
	thisSB->stackDepth = thisSB->stackDepth + 1;
	return;
}

void pop(struct stackBase * thisSB)
{
	struct stackNode * temp = thisSB->firstNode;
	if (temp == NULL)
	{
		return;
	}
	else
	{
		thisSB->firstNode = temp->nextNode;
		thisSB->stackDepth = thisSB->stackDepth - 1;
		free(temp->pointArray);
		free(temp);
		return;
	}
}

int peekDepth(struct stackBase * thisSB)
{
	return (thisSB->stackDepth);
}

void clearStack(struct stackBase *thisSB)
{
	while (thisSB->firstNode != NULL)
	{
		pop(thisSB);
	}
	return;
}
