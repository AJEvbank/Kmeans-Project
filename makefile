all:
	mpicc -Wall CommandLineArgs.c DEBUG.c GetKCentroids.c InitializeKM.c Source.c -o source -lm
