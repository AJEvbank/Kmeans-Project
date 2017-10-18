Team Member Names: Austin Eubank, Frank Conlon, Zohreh Safari

Installation: Use the make file to compile the program.

Running: To run type mpirun -np ## ./source .
## - Can be either 1, 2, or 4 processes.

Command Line Arguments:
1: dim, default dim = 128
2: ndata, default ndata = 10000 (should be divisible by 4)
3: maximum size of data point values, default max_double = 1000.00
4: number of clusters, default k = 128 *
5: seed multiplier (any positive non-zero integer), default seedMult = 1
   //This is a multiplier for the seed array { 2, 4, 8, 16 } and it is the
   //seed for the query point.

Note: Command line arguments cannot be skipped.

Frank Ibem's clusters.py can be used. The program will create the file
data.txt iff dim = 2. If this file has been created, then run clusters.py
(included) in the program directory to view a png of the kmeans and its
cluster assignments.

Output: The program will output the parameters used to construct the kmeans
structure and the query point. It will also output the result(s) of the
kmeans search algorithm and the results of the brute force search algorithm.
The brute force algorithm will select only one point. If this point is one
of the search results, then the program will print, "THE RESULT IS CORRECT".
If the point is not one of the search results, then the program will print,
"THE RESULT IS NOT CORRECT". The program will also output the number of points
checked by the kmeans search algorithm.

* - This value must be divisible by four in the parallelized version, but in
this version it is not necessary for it to be divisible by four.
