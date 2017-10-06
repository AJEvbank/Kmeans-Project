all:
	mpicc -Wall CommandLineArgs.c DEBUG.c GetKCentroids.c InitializeKM.c Source.c ClusterizeKM.c -o source -lm
