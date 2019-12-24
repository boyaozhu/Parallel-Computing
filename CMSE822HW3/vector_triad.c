//-----------------------------------------------------------------------------
// vector_triad.c
//
//      Implements performance testing using floating point
//       multiplication & addition.
//
// Authored by Michael Crockatt
// Modified by H. Metin Aktulga
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <sys/time.h>
#include <time.h>

// Retrieves current time.
void get_walltime(double* wcTime) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);
}

void dummy( double * A, double * C, double *D ) {  A[0] = D[0];  }

//-----------------------------------------------------------------------------
// vec_timeit
//      Performs the vector triad operation A = B + C*D element-wise on vectors
//      of length N 'repetitions' times, recording the runtime required.
//
// INPUT:
//      repetitions     - Number of vector triads to perform.
//      N               - Use vectors of length N.
//
// RETURN:
//      Returns the time elapsed (in seconds) as given by the function
//      get_walltime.
//-----------------------------------------------------------------------------
double vec_timeit (int repetitions, int N) {

  int i, rep;
  double start_time, stop_time;

  double * A = (double*) malloc( N * sizeof(double) );
  double * C = (double*) malloc( N * sizeof(double) );
  double * D = (double*) malloc( N * sizeof(double) );

  // Set values of vector.
  for (i = 0; i < N; ++i) {
    A[i] = 0.0;
    C[i] = rand()%100;
    D[i] = rand()%100;
  }

  get_walltime( &start_time );

  for (rep = 0; rep < repetitions; ++rep) {
    // Compute A += C*D.
    for (i = 0; i < N; ++i) {
      A[i] += C[i]*D[i];
    }
    dummy(A,C,D);
  }

  get_walltime( &stop_time );

  free(A);
  free(C);
  free(D);

  // Return runtime.
  return stop_time - start_time;
}


//-----------------------------------------------------------------------------
// int main
//-----------------------------------------------------------------------------
int main () {

  int i, N, R, ops;
  double elapsed_time;

  // --- Code to compute performance for multiple dimension values. ------ //
  int num = 100;

  srand(time(NULL));

  for ( i = 10; i <= num; ++i ) {
    N = (int)(pow( 10, 7.0/num * i )) / 16 * 16;
    if (N < 1) continue;

    R = 1e8 / N;

    ops = 2 * N * R;

    elapsed_time = vec_timeit( R, N );

    printf( "%10d  %10d  %.3f  %.3f\n",
	    N, R, elapsed_time, (double) ops / elapsed_time / 1e9 );
  }
}
