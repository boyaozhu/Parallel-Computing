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
    
    double *input_array;   // the input array
    double *bucketlist;    //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
    double *local_array;    //This array will contain the elements in each process
    double *total_bucket;
    
    int n, nbar;
    int p;
    int i;
    int my_rank;

    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *bin_elements;       //it will keep track of how many elements have been included in the pth bin

    MPI_Comm_size(MPI_COMM_WORLD, &p);            //processor count
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    //processor rank
    
    n = atoi(argv[1]);    //total elements to be sorted
    nbar = (int)(n/p);      //total elements to be sorted in each process
    n = nbar * p;

    
    double tstart, tend;
    
    tstart = MPI_Wtime();
    input_array = malloc(nbar*sizeof(double));
    //bucketlist = malloc(nbar*sizeof(double));
    bucketlist = calloc(nbar,sizeof(double));
    scounts = malloc(p*sizeof(int));
    dspls = malloc(p*sizeof(int));
    bin_elements = malloc(p*sizeof(int));

    
    srand(time(0)+my_rank);
    double tmp
    for(i = 0 ; i < nbar ; i++){
        tmp = ((double)rand()/(double)(RAND_MAX+1.0));
        input_array[i] = tmp*tmp;
    }
    

    for(i = 0 ; i < p ; i++){
        scounts[i] = 0 ;
    }

    //counting the elements in each processor
    for(i = 0 ; i < nbar ; i++){
        scounts[(int)(input_array[i]/(1.0/p))]++;
    }

    for(i = 0 ; i<p ; i++){
        bin_elements[i] = scounts[i];
    }
    
    
    dspls[0] = 0;
    for(i = 0 ; i< p-1 ;i++){
        dspls[i+1] = dspls[i] + scounts[i];
    }
    
    
    int bin;
    int pos;
    for(i = 0 ; i < nbar ; i++){
        bin = (int)(input_array[i]/(1.0/p));
        pos = dspls[bin] + scounts[bin] - bin_elements[bin];
        bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }
   
    // each process sorts its buckets
    //--------------------------------------------------------//
    /*
    if (my_rank == 0){
        for (i = 0; i < nbar; i++){
            printf("%.2f  ",bucketlist[i]);
        }
        printf("\n");
    }
    else{
        for (i = 0; i < nbar; i++)
            printf("%.2f  ",bucketlist[i]);
        printf("\n");
    }
    */
    int *local_scounts = (int *)malloc(p*sizeof(int));
    
    MPI_Allreduce(scounts, local_scounts, p, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    double *local_bucket = (double *)malloc(local_scounts[my_rank]*sizeof(double));
    
    int proc_id;
    int *scounts_temp;
    int *dsplss;
    for (proc_id = 0; proc_id < p ; proc_id++){
        scounts_temp = malloc(p*sizeof(int));
        MPI_Gather(&scounts[proc_id],1,MPI_INT,scounts_temp,1,MPI_INT,proc_id,MPI_COMM_WORLD);
        MPI_Bcast(scounts_temp,p,MPI_INT,proc_id,MPI_COMM_WORLD);
        if (my_rank == proc_id){
            dsplss = malloc(p*sizeof(int));
            dsplss[0] = 0;
            for(i = 0 ; i< p-1 ;i++){
                dsplss[i+1] = dsplss[i] + scounts_temp[i];
            }
        }
        
        
        MPI_Gatherv(bucketlist ,scounts_temp[my_rank],MPI_DOUBLE,local_bucket,scounts_temp,dsplss,MPI_DOUBLE,proc_id,MPI_COMM_WORLD);
        
        if (my_rank == proc_id){
            free(dsplss);
        }
    
        /*
        nbar = nbar - scounts[proc_id];
        //temp = malloc(nbar*sizeof(double));
        temp = &bucketlist[scounts[proc_id]];
        double *bucketlist_temp = malloc(nbar*sizeof(double));
        bucketlist = temp;
        //bucketlist = bucketlist + scounts[proc_id];
        */
        bucketlist = bucketlist + scounts[proc_id];
        free(scounts_temp);
    }
    /*
    if (my_rank == 0){
        for (i = 0; i < local_scounts[0]; i++){
            printf("%.2f  ",local_bucket[i]);
        }
        printf("\n");
    }
    else{
        for (i = 0; i < local_scounts[1]; i++)
            printf("%.2f  ",local_bucket[i]);
        printf("\n");
    }
    */
    
    
    qsort_dbls(local_bucket,local_scounts[my_rank]);
    
    int *dsplsss;
    if (my_rank == 0){
        //int *dsplss;
        total_bucket = (double *)malloc(n*sizeof(double));
        dsplsss = malloc(p*sizeof(int));
        dsplsss[0] = 0;
        for(i = 0 ; i< p-1 ;i++){
            dsplsss[i+1] = dsplsss[i] + local_scounts[i];
        }
        
    }
    
    
    
    MPI_Gatherv(local_bucket,local_scounts[my_rank],MPI_DOUBLE,total_bucket,local_scounts,dsplsss,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    
    
    tend = MPI_Wtime();
    
    
    double *total_array;
    if (my_rank==0){
        total_array = malloc(n*sizeof(double));
    }
    
    MPI_Gather(input_array,nbar,MPI_DOUBLE,total_array,nbar,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    
    if (my_rank == 0){
        printf("Bining took %.5f seconds\n", tend - tstart);
        qsort_dbls(total_array,n);
        for ( i = 0; i < n; i++) {
            if (total_array[i] != total_bucket[i]) {
                printf("Incorrect result!\n");
                break;
            }
            if (i == n-1) {
                printf("Correct result!\n");
            }
        }
        free(total_array);
    }
    
    if (my_rank == 0){
        free(total_bucket);
        free(dsplsss);

    }
    
    free(dspls);
    free(input_array);
    //free(bucketlist);
    free(bin_elements);
    free(scounts);
    free(local_bucket);
    free(local_scounts);
    //if (my_rank==0) free(total_array);
    MPI_Finalize();
    
    
    /*
     Each bucket can now be sorted individually, either using a different sorting algorithm (i.e., qsort_dbls())
     as you are asked to do, or by recursively applying the bucket sorting algorithm, but we leave this part for
     the parallel code.
     Note that the scount array contains the number of elements each process will receive from the root, dspls
     contains the relative offsets where elements from different processes will be stored, and bucketlist array
     contains the elements stored in the order of the processors.
     */
}
