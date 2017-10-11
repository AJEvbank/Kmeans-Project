#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ROOT 0
#define NUM_ELEMENTS 

int main( int argc, char *argv[] )
{
  //Initialize MPI
  MPI_Init( &argc, &argv );

  //rank: this processe's rank
  //size: the number of processes being used
  //rand_num: the randomly generated number on this process
  //bcast_array: the send buffer for the root process and receive buffer for all others
  int rank, i, *bcast_array;

  //get the rank and size of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  //initialize the recieve buffer on the root  
  bcast_array = (int*)malloc(sizeof(int)*NUM_ELEMENTS);
  
  if(rank == ROOT){
    for(i = 0; i < NUM_ELEMENTS; i++){
      bcast_array[i] = rand();
    }
  }
  else{
    for(i = 0; i < NUM_ELEMENTS; i++){
      bcast_array[i] = 0;
    }
  }

  //have each process print out the numbers it generated
  //This isn't synchronized for this example and may not work well
  printf("Before broadcast the numbers on process %d are:\n", rank);
  for(i = 0; i < NUM_ELEMENTS; i++){
    printf("  %d\n", bcast_array[i]);
  }

  //bcasts data from one process to all processes
  //the i*num_elements index will be the start of the data from the i'th process
  MPI_Bcast(bcast_array, //the send buffer that contains the data to send. This is a pointer
             NUM_ELEMENTS, //The number of items in the buffer to send
             MPI_INT, //The MPI type of the items in the buffer
             ROOT, //the rank of the process to send the bcast
             MPI_COMM_WORLD); //The MPI_Comm to use

  //wait for processes to finish printing their numbers
  MPI_Barrier(MPI_COMM_WORLD);
  //iterate through the recieved elements on each process
  printf("After broadcast the numbers on process %d are:\n", rank);
  for(i = 0; i < NUM_ELEMENTS; i++){
    printf("  %d\n", bcast_array[i]);
  }
  
  free(bcast_array);

  MPI_Finalize();
  return 0;
}
