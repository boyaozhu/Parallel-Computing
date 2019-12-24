
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#include <sys/time.h> // for clock_gettime()

#include "naiveMultiplication.h"
#include "optMultiplication.h"

//#define INT_MAX 100000
int NVEC=16;

void initMatrix(double *matrix, int N, int M)
{
    srand(time(NULL));

    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            matrix[i*M+j] = rand()*1.0/INT_MAX;
            //            printf("%.2lf ", matrix[i*M+j]);
        }
        //        printf("\n");
    }
}

void clear_output_matrix(double* matrix,int N){
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < NVEC; ++j ) {
            matrix[i*NVEC + j] = 0.0;
        }
    }
}

void initVector(double *matrix, int M)
{
    srand(time(NULL));

    for ( int i = 0; i < M; ++i ) {
        for ( int j = 0; j < NVEC; ++j ) {
            matrix[i*NVEC+j] = rand()*1.0/INT_MAX;
        }
    }

}




double getElapsed(struct timeval *start, struct timeval *end)
{
    double secs_used=(end->tv_sec - start->tv_sec); //avoid overflow by subtracting first
    double micros_used= ((secs_used*1000000) + end->tv_usec) - (start->tv_usec);

    return micros_used/1000000.0;
}

int checkDifference(double *matrix1, double *matrix2, int N)
{
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < NVEC; ++j ) {
            if (matrix1[i*NVEC+j] != matrix2[i*NVEC+j]) {
                fprintf(stderr, "matrices differ at position [%d]: ", j);
                fprintf(stderr, "matrix1[%d] = %.5f, ",  j, matrix1[j]);
                fprintf(stderr, "matrix2[%d] = %.5f\n",  j, matrix2[j]);
                return 0;
            }
        }
    }


    return 1;
}

int main(int argc, char* argv[]){
    // Variable declarations

    //  double* naiveinput, *naiveoutput, *optinput, *optoutput;
    int  N, M, B;
    double* input_matrix,*input_vector,*naiveoutput,*optoutput;

    struct timeval start, end;
    double normal_time , blocked_time;
    FILE *fin;

    int testing = 0;
    int num_test = 0;
    int curr_test = 0;
    int multiply_times = 5;

    double total_naive_time = 0.0;
    double total_opt_time = 0.0;



    if (argc == 3 || argc == 5) {

        if (!strcmp(argv[1],"test")) {
            testing = 1;
            fin = fopen(argv[2], "r");
            fscanf(fin, " %d", &num_test);
        } else if (!strcmp(argv[1],"perf")) {
            testing = 0;
        } else {
            fprintf(stderr, "\nERROR: Unknown option %s\n", argv[1]);
        }

    } else {
        fprintf(stderr, "\nERROR: This program can be executed in the \
                \"test\" or \"perf\" modes:\n\n");
        fprintf(stderr,"           multiplication test inputfile\n");
        fprintf(stderr,"           multiplication perf N M B\n\n");
        fprintf(stderr,"           inputfile format: starts with a line indicating \
                the number of tests,\n");
        fprintf(stderr,"                             followed the tests specified with \
                N M B values.\n\n");
        fprintf(stderr,"           N: number of matrix rows\n");
        fprintf(stderr,"           M: number of matrix columns\n");
        fprintf(stderr,"           B: blocking factor.\n");
        return 0;
    }

    // If doing testing run, print header for output to CSV file.
    if (testing) printf( "N,M,B,naive_time,opt_time,opt_correct\n" );


    // main loop: executed only once for performance testing to ensure cold start.
    //            However, it is okay to do several accuracy checks (test runs)
    //            at once in a single execution.



    do {

        if (testing) {
            curr_test++;
            // read from the input file
            fscanf(fin, " %d %d %d", &N, &M, &B);
        } else {
            // otherwise get the command line arguments
            N = strtol(argv[2], NULL, 10);
            M = strtol(argv[3], NULL, 10);
            B = strtol(argv[4], NULL, 10);
            //            curtest++;
        }

        if (testing) {
            printf( "%d,%d,%d,", N, M, B );
        } else {
            printf("ninps = %d, curtest = %d\n", num_test, curr_test);
            printf("N = %d, M = %d, B = %d\n", N, M, B);
        }

        //Allocations
        input_matrix = (double *)malloc(N*M*sizeof(double));
        input_vector = (double *)malloc(M*NVEC*sizeof(double));
        naiveoutput = (double *)malloc(N*NVEC*sizeof(double));
        optoutput = (double *)malloc(N*NVEC*sizeof(double));


        //initializing the matrix and vector
        initMatrix(input_matrix,N,M);
        initVector(input_vector,M);

        //run the both versions multiple times and take their average
        for(int t = 0; t < multiply_times ; t++){
            //naiveMultiplication
            clear_output_matrix(naiveoutput,N);
            gettimeofday(&start,NULL);
            naiveMultiplication(naiveoutput,input_matrix,input_vector,N,M);
            gettimeofday(&end,NULL);
            normal_time = getElapsed(&start,&end);


            //blocked multiplication
            clear_output_matrix(optoutput,N);
            gettimeofday(&start,NULL);
            optMultiplication(optoutput,input_matrix,input_vector,N,M,B);
            gettimeofday(&end,NULL);
            blocked_time = getElapsed(&start,&end);


            // checkDifference(naiveoutput,optoutput,N);

            total_naive_time += normal_time;
            total_opt_time += blocked_time;
        }

        // How do results compare?
        if (testing) {
            printf( "%.4e,%.4e,%d\n", total_naive_time/multiply_times, total_opt_time/multiply_times,
                    checkDifference(naiveoutput, optoutput, N) );
        } else {
            // checkDifference(naiveoutput,optoutput,N);
            printf("Naive time: %.3f sec\n", total_naive_time/multiply_times);
            printf("Opt. time:  %.3f sec\n", total_opt_time/multiply_times);
            printf("Achieved speedup= %.2f\n\n", total_naive_time/total_opt_time);
        }
        
        total_naive_time = 0.0;
        total_opt_time = 0.0;

        free(input_matrix);
        free(input_vector);
        free(naiveoutput);
        free(optoutput);
    } while (curr_test<num_test);






    return 0;
}
