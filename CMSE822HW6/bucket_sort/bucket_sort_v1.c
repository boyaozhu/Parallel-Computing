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
    
    
    MPI_Init(&argc, &argv);
    double *input_array;    //the input array
    double *bucketlist;     //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
	double *local_bucket;   //contain elements of each processors
    double *local_array;    //This array will contain the elements in each process
	double t_start, t_gen, t_bin, t_dis, t_lsrt, t_gthr, t_tot, t_temp;    //the sys time
    int n, i, p, my_rank;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *bin_elements;       //it will keep track of how many elements have been included in the pth bin

	MPI_Comm_size(MPI_COMM_WORLD, &p);			//processor count
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);	//processor rank

	t_start = MPI_Wtime();

	if (my_rank == 0)
    {
        n = atoi(argv[1]);
        
        printf("%d   %d",n,p);
        printf("\n\n");
        
		input_array = malloc(n * sizeof(double));
		bucketlist = malloc(n * sizeof(double));
		scounts = malloc(p * sizeof(int));
		dspls = malloc(p * sizeof(int));
		bin_elements = malloc(p * sizeof(int));

		srand(time(0)+my_rank);
		// generate

		t_start = MPI_Wtime();
		for (i = 0; i < n; i++) {
			input_array[i] = ((double)rand() / (double)(RAND_MAX));
		}
		t_gen = MPI_Wtime() - t_start;
        
		// bin
		// initialize counts in each processors
		t_temp = MPI_Wtime();

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
    
	if (my_rank == 0) {
		t_bin = MPI_Wtime() - t_temp;
	}
    
	// distribute
	t_temp = MPI_Wtime();
	if (my_rank != 0) {
		scounts = malloc(p * sizeof(int));
	}
	MPI_Bcast(scounts, p, MPI_INT, 0, MPI_COMM_WORLD);
	// send bucket to other processors
	local_bucket = malloc(scounts[my_rank] * sizeof(double));
    
	MPI_Scatterv(bucketlist, scounts, dspls, MPI_DOUBLE, local_bucket, scounts[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
	t_dis = MPI_Wtime() - t_temp;

	// local sort
	t_temp = MPI_Wtime();
	qsort_dbls(local_bucket, scounts[my_rank]);
	t_lsrt = MPI_Wtime() - t_temp;
  
	// collect result
	t_temp = MPI_Wtime();

	MPI_Gatherv(local_bucket, scounts[my_rank], MPI_DOUBLE, bucketlist, scounts, dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
	t_gthr = MPI_Wtime() - t_temp;

	// total time
	t_tot = MPI_Wtime() - t_start;

	if (my_rank == 0) {
		// Generate, Bin, Distribute, Local Sort, Gather, Total
		printf("Generate = %.6f\n Bin = %.6f\n Distribute = %.6f\n Local Sort = %.6f\n Gather = %.6f\n Total = %.6f\n",
			t_gen, t_bin, t_dis, t_lsrt, t_gthr, t_tot);
        // verification
        
        for ( i = 0; i < n-1; i++) {
            if (bucketlist[i] > bucketlist[i+1]) {
                printf("Incorrect result!\n");
                break;
            }
            if (i == n-2) {
                printf("Correct result!\n");
            }
        }
        
        free(input_array);
        free(bucketlist);
        free(dspls);
        free(bin_elements);
        printf("\n\n");
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
