#include "optMultiplication.h"
#include <stdio.h>
void optMultiplication(double* output, double* input_matrix, double* input_vector,int N,int M,int B){
    //your code here
    int i, j, k, ii, jj;
    int enR = B * (N/B);   //Amount that fits evenly into row
    int enC = B * (M/B);   //Amount that fits evenly into column
    
    int reR = N - enR;
    int reC = M - enC;
    
    
    for (ii = 0; ii < enR; ii += B){ //go through each block through row, increment by B
        for (jj = 0; jj < enC; jj += B){ //go through each block through column, increment by B
            for (i = ii; i < ii + B; i++){ //go through each element to the right for each block
                for (j = jj; j < jj + B; j++){ //go through each element down for each block
                    for (k = 0; k < 16; k++){ //go through each element to the right
                        output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
        
                    }
                }
            }
        }
    }
    
    //the following code implements when blocks don't evenly it the matrix
    if ((reR != 0)&&(reC == 0)){
        for (i = enR; i < N; i++){
            for (j = 0; j < M; j++){
                for (k = 0; k < 16; k++){
                    output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
                }
            }
        }
    }
        
    if ((reR == 0)&&(reC != 0)){
        for (i = 0; i < N; i++){
            for (j = enC; j < M; j++){
                for (k = 0; k < 16; k++){
                    output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
                }
            }
        }
    }
    
    if ((reR != 0)||(reC != 0)){
        for (i = enR; i < N; i++){
            for (j = 0; j < enC; j++){
                for (k = 0; k < 16; k++){
                    output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
                }
            }
        }
        for (i = 0; i < enR; i++){
            for (j = enC; j < M; j++){
                for (k = 0; k < 16; k++){
                    output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
                }
            }
        }
        for (i = enR; i < N; i++){
            for (j = enC; j < M; j++){
                for (k = 0; k < 16; k++){
                    output[i*16+k] += input_matrix[i*M+j]*input_vector[j*16+k];
                }
            }
        }
        
    }
    
}

