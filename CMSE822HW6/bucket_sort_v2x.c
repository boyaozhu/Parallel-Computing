#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


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

    double *input_array;   // the input array
    double *bucketlist;    //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
    double *local_bucketlist;
    double *local_array;    //This array will contain the elements in each process

    double tmp, t_start, t_gen, t_bin, t_dis, t_lsrt, t_gthr, t_tot, t_temp;

    int n, p, i, my_rank, local_size, local_n;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *recv_scounts;      //scounts for receiving from other processes
    int *local_scounts;
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *recv_dspls;        //dspls for receiving from other processes
    int *local_dspls;
    int *bin_elements;       //it will keep track of how many elements have been included in the pth bin

    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    n = atoi(argv[1]);
    local_n = n / p;
    if (my_rank == 0) {
        local_n += n % p;
    }
        

    input_array = malloc(local_n * sizeof(double));
    bucketlist = malloc(n * sizeof(double));
    local_bucketlist = malloc(local_n * sizeof(double));
    scounts = malloc(p * sizeof(int));
    recv_scounts = malloc(p * sizeof(int));
    local_scounts = malloc(p * sizeof(int));
    dspls = malloc(p * sizeof(int));
    recv_dspls = malloc(p * sizeof(int));
    local_dspls = malloc(p * sizeof(int));
    bin_elements = malloc(p * sizeof(int));

    // seed with time
    srand(time(0)+my_rank);

    t_start = MPI_Wtime();

    // generate
    for(i = 0 ; i < local_n ; i++){
        tmp = ((double) rand()/RAND_MAX);
        input_array[i] = tmp;
	// input_array[i] = ((double) rand()/RAND_MAX);
    }
    t_gen = MPI_Wtime() - t_start;

    // bin
    t_temp = MPI_Wtime();

    for(i = 0 ; i < p ; i++){
        local_scounts[i] = 0 ;
    }

    for(i = 0 ; i < local_n ; i++){
        local_scounts[(int)(input_array[i]/(1.0/p))]++;
    }

    for(i = 0 ; i<p ; i++){
        bin_elements[i] = local_scounts[i];
    }

    local_dspls[0] = 0;
    for(i = 0 ; i < p-1 ;i++){
        local_dspls[i+1] = local_dspls[i] + local_scounts[i];
    }

    int bin;
    int pos;
    for(i = 0 ; i < local_n ; i++){
        bin = (int)(input_array[i]/(1.0/p));
        pos = local_dspls[bin] + local_scounts[bin] - bin_elements[bin];
        local_bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }

    t_bin = MPI_Wtime() - t_temp;

    // distribute
    t_temp = MPI_Wtime();

    // figure out receving counts and offsets
    MPI_Alltoall(local_scounts, 1, MPI_INT, recv_scounts, 1, MPI_INT, MPI_COMM_WORLD);

    recv_dspls[0] = 0;
    for (i = 1; i < p; i++){
        recv_dspls[i] = recv_dspls[i-1] + recv_scounts[i-1];
    }

    local_size = recv_dspls[p-1] + recv_scounts[p-1];

    local_array = malloc(local_size*sizeof(double));
    MPI_Alltoallv(local_bucketlist, local_scounts, local_dspls, MPI_DOUBLE, 
                  local_array, recv_scounts, recv_dspls, MPI_DOUBLE, MPI_COMM_WORLD);

    t_dis = MPI_Wtime() - t_temp;

    // local sort
    t_temp = MPI_Wtime();

    qsort_dbls(local_array, local_size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    t_lsrt = MPI_Wtime() - t_temp;

    // gather
    t_temp = MPI_Wtime();

    MPI_Reduce(local_scounts, scounts, p, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0){
        dspls[0] = 0;
        for (i = 1; i < p; i++){
            dspls[i] = dspls[i-1] + scounts[i-1];
        }
    }

    MPI_Gatherv(local_array, local_size, MPI_DOUBLE, 
                bucketlist, scounts, dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    t_gthr = MPI_Wtime() - t_temp;

    // total time
    t_tot = MPI_Wtime() - t_start;


    // validating sorted array
    if (my_rank == 0){
        for (i = 0; i < n - 1; i++)
            if (bucketlist[i] > bucketlist[i+1])
                printf("Incorrect sorting found at position %d/%d, %.3f > %.3f\n", 
                       i+1, n, bucketlist[i], bucketlist[i+1]);
        // Generate, Bin, Distribute, Local Sort, Gather, Total
        printf("Number_of_processors = %d\nGenerate = %.6f\nBin = %.6f\nDistribute = %.6f\nLocal_sort = %.6f\nGather = %.6f\nTotal = %.6f\n",
               p, t_gen, t_bin, t_dis, t_lsrt, t_gthr, t_tot);
    }

    MPI_Finalize();
}
