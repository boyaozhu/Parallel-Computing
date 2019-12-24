/* 
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory 
 * and reimplementation by Scott B. Baden, UCSD
 * 
 * Modified and  restructured by Didem Unat, Koc University
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <cuda.h>
using namespace std;


// Utilities
// 

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
int blocksize = 32; 

double getTime()
{
    struct timeval TV;
    struct timezone TZ;

    const int RC = gettimeofday(&TV, &TZ);
    if(RC == -1) {
            cerr << "ERROR: Bad call to gettimeofday" << endl;
            return(-1);
    }

    return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );

}  // end getTime()

// Allocate a 2D array
double *alloc2D(int m,int n){
   double *E;
   int nx=n, ny=m;
   E = (double*)malloc(sizeof(double)*nx*ny);
   assert(E);
   return(E);
}
    
// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
 double stats(double *E, int m, int n, double *_mx){
     double mx = -1;
     double l2norm = 0;
     int i, j;
     for (j=1; j<=m; j++)
       for (i=1; i<=n; i++) {
	   l2norm += E[j*(n+2)+i]*E[j*(n+2)+i];
	   if (E[j*(n+2)+i] > mx)
	       mx = E[j*(n+2)+i];
      }
     *_mx = mx;
     l2norm /= (double) ((m)*(n));
     l2norm = sqrt(l2norm);
     return l2norm;
 }

// External functions
extern "C" {
    void splot(double *E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double& T, int& n, int& px, int& py, int& plot_freq, int& no_comm, int&num_threads);

__global__ void setBoundaryLeft(double* E_prev, const int m, const int n)
{
    int j = blockIdx.x*blockDim.x+threadIdx.x+1; 
    if (j<=m) E_prev[j*(n+2)+0] = E_prev[j*(n+2)+2];
}
__global__ void setBoundaryRight(double* E_prev, const int m, const int n)
{
    int j = blockIdx.x*blockDim.x+threadIdx.x+1; 
    if (j<=m) E_prev[j*(n+2)+n+1] = E_prev[j*(n+2)+n-1];
}
__global__ void setBoundaryUp(double* E_prev, const int m, const int n)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x+1; 
    if (i<=n) E_prev[i] = E_prev[2*(n+2)+i];
}
__global__ void setBoundaryDown(double* E_prev, const int m, const int n)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x+1; 
    if (i<=n) E_prev[(m+1)*(n+2)+i] = E_prev[(m-1)*(n+2)+i];
}

__global__ void solvePDE(double *E, double *E_prev, const double alpha, const int m, const int n)
{
    int j = blockIdx.y*blockDim.y+threadIdx.y+1; 
    int i = blockIdx.x*blockDim.x+threadIdx.x+1;
    
    if ((j<=m) && (i<=n))
    {
        E[j*(n+2)+i] = E_prev[j*(n+2)+i]+alpha*(E_prev[j*(n+2)+i+1]+E_prev[j*(n+2)+i-1]-4*E_prev[j*(n+2)+i]+E_prev[(j+1)*(n+2)+i]+E_prev[(j-1)*(n+2)+i]); 
    }
}

__global__ void solveODE_E(double *E, double *R,  const double kk,
	       const double dt, const double a, const int m, const int n)
{
    int j = blockIdx.y*blockDim.y+threadIdx.y+1; 
    int i = blockIdx.x*blockDim.x+threadIdx.x+1;
    
    if ((j<=m) && (i<=n))
    {
        E[j*(n+2)+i] = E[j*(n+2)+i] -dt*(kk* E[j*(n+2)+i]*(E[j*(n+2)+i] - a)*(E[j*(n+2)+i]-1)+ E[j*(n+2)+i] *R[j*(n+2)+i]);
    }
}

__global__ void solveODE_R(double *E, double *R,  const double kk,
	       const double dt, const double epsilon, const double M1,
           const double  M2, const double b, const int m, const int n)
{
    int j = blockIdx.y*blockDim.y+threadIdx.y+1; 
    int i = blockIdx.x*blockDim.x+threadIdx.x+1;
    
    if ((j<=m) && (i<=n))
    {
        R[j*(n+2)+i] = R[j*(n+2)+i] + dt*(epsilon+M1* R[j*(n+2)+i]/( E[j*(n+2)+i]+M2))*(-R[j*(n+2)+i]-kk* E[j*(n+2)+i]*(E[j*(n+2)+i]-b-1));
    }
}


void simulate (double* E,  double* E_prev,double* R,
	       const double alpha, const int n, const int m, const double kk,
	       const double dt, const double a, const double epsilon,
	       const double M1,const double  M2, const double b)
{

    /* 
     * Copy data from boundary of the computational box 
     * to the padding region, set up for differencing
     * on the boundary of the computational box
     * Using mirror boundaries
     */
      
    setBoundaryLeft<<<ceil(m/64.0),64>>>(E_prev, m, n); 
    setBoundaryRight<<<ceil(m/64.0),64>>>(E_prev, m, n);
    setBoundaryUp<<<ceil(n/64.0),64>>>(E_prev, m, n);
    setBoundaryDown<<<ceil(n/64.0),64>>>(E_prev, m, n);
    
    
    // Solve for the excitation, the PDE
    dim3 DimGrid(ceil((double)n/blocksize), ceil((double)m/blocksize), 1); 
    dim3 DimBlock(blocksize, blocksize, 1); 
    solvePDE<<<DimGrid, DimBlock>>>(E, E_prev, alpha, m, n); 
    
    /* 
     * Solve the ODE, advancing excitation and recovery to the
     *     next timtestep
     */
    solveODE_E<<<DimGrid, DimBlock>>>(E, R, kk, dt, a, m, n); 
    solveODE_R<<<DimGrid, DimBlock>>>(E, R, kk, dt, epsilon, M1, M2, b, m, n); 
    
}

// Main program
int main (int argc, char** argv)
{
  /*
   *  Solution arrays
   *   E is the "Excitation" variable, a voltage
   *   R is the "Recovery" variable
   *   E_prev is the Excitation variable for the previous timestep,
   *      and is used in time integration
   */
  double *E, *R, *E_prev;
  
  // Various constants - these definitions shouldn't change
  const double a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;
  
  double T=1000.0;
  int m=200,n=200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads=1; 

  cmdLine( argc, argv, T, n,px, py, plot_freq, no_comm, num_threads);
  m = n;  
  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the 
  // boundaries of the computation box
  E = alloc2D(m+2,n+2);
  E_prev = alloc2D(m+2,n+2);
  R = alloc2D(m+2,n+2);
  
  int i,j;
  // Initialization
  for (j=1; j<=m; j++)
    for (i=1; i<=n; i++)
      E_prev[j*(n+2)+i] = R[j*(n+2)+i] = 0;
  
  for (j=1; j<=m; j++)
    for (i=n/2+1; i<=n; i++)
      E_prev[j*(n+2)+i] = 1.0;
  
  for (j=m/2+1; j<=m; j++)
    for (i=1; i<=n; i++)
      R[j*(n+2)+i] = 1.0;
  
  double dx = 1.0/n;

  // For time integration, these values shouldn't change 
  double rp= kk*(b+1)*(b+1)/4;
  double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
  double dtr=1/(epsilon+((M1/M2)*rp));
  double dt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
  double alpha = d*dt/(dx*dx);

  cout << "Grid Size       : " << n << endl; 
  cout << "Duration of Sim : " << T << endl; 
  cout << "Time step dt    : " << dt << endl; 
  cout << "Process geometry: " << px << " x " << py << endl;
  if (no_comm)
    cout << "Communication   : DISABLED" << endl;
  
  cout << endl;
  
  //cuda: transfer data to GPU
  double *E_cuda, *R_cuda, *E_prev_cuda;
  int size = (m+2)*(n+2)*sizeof(double);
  cudaMalloc((void **) &E_cuda, size);
  cudaMemcpy(E_cuda, E, size, cudaMemcpyHostToDevice); 
  cudaMalloc((void **) &R_cuda, size);
  cudaMemcpy(R_cuda, R, size, cudaMemcpyHostToDevice); 
  cudaMalloc((void **) &E_prev_cuda, size);
  cudaMemcpy(E_prev_cuda, E_prev, size, cudaMemcpyHostToDevice); 
  
  // Start the timer
  double t0 = getTime();
  
 
  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;
  // Integer timestep number
  int niter=0;
  
  while (t<T) {
    
    t += dt;
    niter++;
 
    simulate(E_cuda, E_prev_cuda, R_cuda, alpha, n, m, kk, dt, a, epsilon, M1, M2, b); 
    
    //swap current E with previous E
    double *tmp = E_cuda; E_cuda = E_prev_cuda; E_prev_cuda = tmp;
    
    if (plot_freq){
      int k = (int)(t/plot_freq);
      if ((t - k * plot_freq) < dt){
        cudaMemcpy(E, E_cuda, size, cudaMemcpyDeviceToHost);  
	    splot(E,t,niter,m+2,n+2);
      }
    }
  }//end of while loop
  
  double time_elapsed = getTime() - t0;

  double Gflops = (double)(niter * (1E-9 * n * n ) * 28.0) / time_elapsed ;
  double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0  ))/time_elapsed;

  cout << "Number of Iterations        : " << niter << endl;
  cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
  cout << "Sustained Gflops Rate       : " << Gflops << endl; 
  cout << "Sustained Bandwidth (GB/sec): " << BW << endl << endl; 

  //cuda: transfer data from GPU
  cudaMemcpy(E_prev, E_prev_cuda, size, cudaMemcpyDeviceToHost); 
  cudaMemcpy(E, E_cuda, size, cudaMemcpyDeviceToHost); 
  cudaMemcpy(R, R_cuda, size, cudaMemcpyDeviceToHost); 
  cudaFree(E_cuda); cudaFree(R_cuda); cudaFree(E_prev_cuda); 
  
  double mx;
  double l2norm = stats(E_prev,m,n,&mx);
  cout << "Max: " << mx <<  " L2norm: "<< l2norm << endl;

  if (plot_freq){
    cout << "\n\nEnter any input to close the program and the plot..." << endl;
    getchar();
  }
  
  free (E);
  free (E_prev);
  free (R);
  
  return 0;
}
