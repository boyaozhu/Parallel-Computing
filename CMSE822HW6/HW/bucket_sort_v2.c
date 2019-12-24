#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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
    
    double t_start, t_end, t1, t2;      
    double *input_array;   // the input array
    double *bucketlist;    //this array will contain input elements in order of the processors
                            //e.g elements of process 0 will be stored first, then elements of process 1, and so on
    double *local_array;    //This array will contain the elements in each process
    
    int n;
    int p;
    int i;
    int my_rank;
    
    int *scounts;           //This array will contain the counts of elements each processor will receive
    int *dspls;             //The relative offsets in bucketlist array where the elements of different processes
                            //will be stored
    int *bin_elements;       //it will keep track of how many elements have been included in the pth bin
    
    int my_ranks, my_n; 
    
    MPI_Init(&argc, &argv); 
    t_start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &p);  // get the total number of processes 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  // get the rank of my process
    
    //******Generate
    t1 = MPI_Wtime(); 
    n = atoi(argv[1]);
    my_n = n/p; 
    if (my_rank + 1 == p) my_n = my_n + n%p;
    
    if (my_rank==0){
        printf("%d   %d",n,p);
        printf("\n\n");
    }
    
    input_array = malloc(my_n*sizeof(double));
    bucketlist = malloc(my_n*sizeof(double));
    scounts = malloc(p*sizeof(int));
    dspls = malloc(p*sizeof(int));
    bin_elements = malloc(p*sizeof(int));
    srand(time(0)+my_rank);
    for(i = 0 ; i < my_n ; i++){
        input_array[i] = ((double) rand()/RAND_MAX+1);
    }
    t2 = MPI_Wtime(); 
    if (my_rank == 0) printf("Generate = %f \n", t2-t1);
    //******Generate
    
    //******Bin
    t1 = MPI_Wtime(); 
    for(i = 0 ; i < p ; i++){
        scounts[i] = 0 ;
    }
    //counting the elements in each processor
    for(i = 0 ; i < my_n ; i++){
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
    for(i = 0 ; i < my_n ; i++){
        bin = (int)(input_array[i]/(1.0/p));
        pos = dspls[bin] + scounts[bin] - bin_elements[bin];
        bucketlist[pos] = input_array[i];
        bin_elements[bin]--;
    }
    t2 = MPI_Wtime(); 
    if (my_rank == 0) printf("Bin = %f \n", t2-t1);
    //******Bin
 
    //******Distribute
    t1 = MPI_Wtime();
    int *rec_scounts;
    int *rec_dspls; 
    int local_len;
    rec_scounts = malloc(p*sizeof(int));
    rec_dspls = malloc(p*sizeof(int));
    MPI_Alltoall(scounts, 1, MPI_INT, rec_scounts, 1, MPI_INT, MPI_COMM_WORLD);
    local_len = 0; 
    for (i = 0; i<p ; i++) {
        local_len += rec_scounts[i]; 
    }
    rec_dspls[0] = 0;
    for (i = 0; i<p-1; i++) {
        rec_dspls[i+1] = rec_dspls[i] + rec_scounts[i]; 
    }
    local_array = malloc(local_len*sizeof(double));
    MPI_Alltoallv(bucketlist, scounts, dspls, MPI_DOUBLE, 
        local_array, rec_scounts, rec_dspls, MPI_DOUBLE, MPI_COMM_WORLD);
    t2 = MPI_Wtime(); 
    if (my_rank == 0) printf("Distribute = %f \n", t2-t1);
    //******Distribute
    
    //******Local sort
    t1 = MPI_Wtime();
    qsort_dbls(local_array, local_len);
    t2 = MPI_Wtime(); 
    if (my_rank == 0) printf("Local sort = %f \n", t2-t1);
    //******Local sort
     
    //******Gather
    t1 = MPI_Wtime();
    double *result; 
    int *gat_scounts; 
    int *gat_dspls; 
    if (my_rank == 0) 
    {
        result = malloc(n*sizeof(double));
        gat_scounts = malloc(p*sizeof(int));
        gat_dspls = malloc(p*sizeof(int));
    }
    MPI_Gather(&local_len, 1, MPI_INT, gat_scounts, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    if (my_rank == 0) 
    {
        gat_dspls[0] = 0;
        for (i = 0; i<p-1; i++) {
            gat_dspls[i+1] = gat_dspls[i] + gat_scounts[i]; 
        }
    }
    MPI_Gatherv(local_array, local_len, MPI_DOUBLE, 
        result, gat_scounts, gat_dspls, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    t2 = MPI_Wtime(); 
    if (my_rank == 0) printf("Gather = %f \n", t2-t1);
    //******Gather
    
    t_end = MPI_Wtime(); 
    if (my_rank == 0) printf("Total = %f \n",t_end - t_start);
    
    //******Verification
    if (my_rank == 0)
    {
        int flag = 1; 
        for (i = 0; i<n-1; i++) 
        {
            if (compare_dbls(result+i, result+i+1)>0)
            {
                flag = 0; 
                break;
            }
        }
        if (flag)
            printf("%s\n", "Correct result!");
        else
            printf("%s\n", "Incorrect");
         
    }
    //******Verification
    
    if (my_rank==0){
        
        printf("\n\n\n");
    }
    
    MPI_Finalize(); 
}
