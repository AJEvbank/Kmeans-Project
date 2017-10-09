#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define ROOT 0
#define NUM_ELEMENTS 2

int main( int argc, char *argv[] )
{
  //Initialize MPI
  MPI_Init( &argc, &argv );

  //rank: this processe's rank
  //size: the number of processes being used
  //rand_num: the randomly generated number on this process
  //reduce_array: the receive buffer for the root process
  int rank, size, i, rand_nums[NUM_ELEMENTS], *reduce_array;

  //get the rank and size of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  
  //initialize the recieve buffer on the root  
  reduce_array = (int*)malloc(sizeof(int)*NUM_ELEMENTS);
  
  //generate the random number
  srand(rank*2+1);
  
  for(i = 0; i < NUM_ELEMENTS; i++){
    rand_nums[i] = rand();
  }

  //have each process print out the numbers it generated
  //This isn't synchronized for this example and won't work well
  printf("The random numbers on process %d are:\n", rank);
  for(i = 0; i < NUM_ELEMENTS; i++){
    printf("  %d\n", rand_nums[i]);
  }

  //reduces data from all processes into an array on one process
  //the 0 element will be the data from the rank 0 process
  //the i*num_elements index will be the start of the data from the i'th process
  MPI_Reduce(rand_nums, //the send buffer that contains the data to send. This is a pointer
             reduce_array, //The buffer to put the reduced items into on the recieving process. This is a pointer
             NUM_ELEMENTS, //The number of items in the buffer to send
             MPI_INT, //The MPI type of the items in the buffer
             MPI_MIN, //The reduce function to use
             ROOT, //the rank of the process to recieve the reduce
             MPI_COMM_WORLD); //The MPI_Comm to use

  //wait for processes to finish printing their numbers
  MPI_Barrier(MPI_COMM_WORLD);
  //iterate through the recieved elements on the root process
  if(rank == ROOT){
    for(i = 0; i < NUM_ELEMENTS; i++){
      printf("The smallest number at index %d is %d\n", i, reduce_array[i]);
    }
  }

  //sum reduce now 
  MPI_Reduce(rand_nums, reduce_array, NUM_ELEMENTS, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);

  if(rank == ROOT){
    for(i = 0; i < NUM_ELEMENTS; i++){
      printf("The sum at index %d is %d\n", i, reduce_array[i]);
    }
  }

  //Max all reduce
  MPI_Allreduce(rand_nums, reduce_array, NUM_ELEMENTS, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  for(i = 0; i < NUM_ELEMENTS; i++){
    printf("The max at index %d on process %d is %d\n", i, rank, reduce_array[i]);
  }

  free(reduce_array);
  MPI_Finalize();
  return 0;
}
