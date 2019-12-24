#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>


// Comparison function used by qsort
int compare_dbls(const void* arg1, const void* arg2)
{
    double a1 = *(double *) arg1;
    double a2 = *(double *) arg2;
    if (a1 < a2) return -1;
    else if (a1 == a2) return 0;
    else return 1;
}
// Sort the array in place
void qsort_dbls(double *array, int array_len)
{
    qsort(array, (size_t)array_len, sizeof(double),
          compare_dbls);
}

void main(int argc, char* argv[]){
    
    MPI_Init(&argc,&argv);
    
    double *input_array;    //the input array
    double *bucketlist;     //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
	double *local_bucket;   //contain elements of each processors
    double *local_array;    //This array will contain the elements in each process
	double tstart, tend;    //the sys time
    int n, i, p, my_rank;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *bin_elements;       //it will keep track of how many elements have been included in the pth bin

	MPI_Comm_size(MPI_COMM_WORLD, &p);			//processor count
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	//processor rank

	//tstart = MPI_Wtime();

	if (my_rank == 0)
    {
        n = atoi(argv[1]);
        tstart = MPI_Wtime();
		input_array = malloc(n * sizeof(double));
		bucketlist = malloc(n * sizeof(double));
		scounts = malloc(p * sizeof(int));
		dspls = malloc(p * sizeof(int));
		bin_elements = malloc(p * sizeof(int));

		// input array
        srand(time(0));
		for (i = 0; i < n; i++) {
			input_array[i] = ((double)rand() / (double)(RAND_MAX+1.0));
		}
        
		// initialize counts in each processors
		for (i = 0; i < p; i++) {
			scounts[i] = 0;
		}

		//count the elements in each processor
		for (i = 0; i < n; i++) {
			scounts[(int)(input_array[i] / (1.0 / p))]++;
		}
          
		for (i = 0; i < p; i++) {
			bin_elements[i] = scounts[i];
		}

		dspls[0] = 0;
		for (i = 0; i < p - 1; i++) {
			dspls[i + 1] = dspls[i] + scounts[i];
		}

		int bin;
		int pos;
		// bucket sort
		for (i = 0; i < n; i++) {
			bin = (int)(input_array[i] / (1.0 / p));
			pos = dspls[bin] + scounts[bin] - bin_elements[bin];
			bucketlist[pos] = input_array[i];
			bin_elements[bin]--;
		}
	}
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (my_rank != 0){
        scounts = malloc(p * sizeof(int));
    }
    
	MPI_Bcast(scounts, p, MPI_INT, 0, MPI_COMM_WORLD);
	// send bucket to other processors
	local_bucket = malloc(scounts[my_rank] * sizeof(double));
    
    
	MPI_Scatterv(bucketlist, scounts, dspls, MPI_DOUBLE, local_bucket, scounts[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
	// do local sort
	qsort_dbls(local_bucket, scounts[my_rank]);
  
	// collect result
	MPI_Gatherv(local_bucket, scounts[my_rank], MPI_DOUBLE, bucketlist, scounts, dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
	tend = MPI_Wtime();

	if (my_rank == 0) {
		printf("Bining took %.5f seconds\n", tend - tstart);
        
        // verification
        qsort_dbls(input_array, n);
        for ( i = 0; i < n; i++) {
            if (input_array[i] != bucketlist[i]) {
                printf("Incorrect result!\n");
                break;
            }
            if (i == n-1) {
                printf("Correct result!\n");
            }
        }
        
        free(input_array);
        free(bucketlist);
        free(dspls);
        free(bin_elements);
    }
    free(local_bucket);
    free(scounts);
    /*
     Each bucket can now be sorted individually, either using a different sorting algorithm (i.e., qsort_dbls())
     as you are asked to do, or by recursively applying the bucket sorting algorithm, but we leave this part for
     the parallel code.
     Note that the scount array contains the number of elements each process will receive from the root, dspls
     contains the relative offsets where elements from different processes will be stored, and bucketlist array
     contains the elements stored in the order of the processors.
     */
    
    
    
	MPI_Finalize();
}
