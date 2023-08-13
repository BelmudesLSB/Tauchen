#include <iostream>
#include <vector> 
#include <cmath>
#include "tauchen.hpp" 

int main(){

    // Parameters needed for the discretization:

    const int N = 9;                    // Number of grid points.
    const int m = 3;                    // Number of max standard deviations that the code allows.
    const double lambda = 0.9;          // Persistence of the process.
    const double sigma_e = 0.1;         // Standard deviation of the process.

    //create memory for the vectors:

    double *y = new double[N];
    double *P = new double[N*N];

    // Call the function and display the results:
    
    tauchen(y,P,m,lambda,sigma_e,N);
    displayV(y,N);
    displayV(P,N*N); 
    return 0;
}