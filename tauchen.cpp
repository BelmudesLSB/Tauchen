#include<cmath> 
#include<iostream>
#include<cstdlib>

using namespace std; 

void displayV(double *v, int M){
    for (int i=0;i<M;i++){
        cout << v[i] << endl;
    }
}

double normalCDF(double x){
    return std::erfc(-x / std::sqrt(2)) / 2;
}

//// This function will create a vector of size N*N. The first N elements will be the first row of the matrix,
// the next N elements will be the second row of the matrix, and so on. element [i * N + j] will represent the
// probability that the state i will move to the state j. Note that for each element in the vector, we have a
// biyection mapping to one pair (i,j). (Element % N) will give us the column, and (Element / N) will give us the row.

int tauchen(double *&y_grid, double *&Prob, const int m, const double lambda, const double sigma_e, const int N){

    // This code follows Tauchen (1986) to discretize an AR(1) process using e-N(0,sigma_e^2).

    // Allocate memory for the vectors:

    y_grid = new double[N];
    Prob = new double[N*N];

    // Construct y_grid:

    double sigma_y = sqrt(pow(sigma_e,2)/(1-pow(lambda,2)));
    y_grid[N-1] = m*sigma_y;
    y_grid[0] = -m*sigma_y;
    double omega = (y_grid[N-1]-y_grid[0])/(N-1);

    for (int i=1; i<N-1; i++){
        y_grid[i] = y_grid[i-1] + omega;
    }

    // Construct Prob:

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            // Following Tauchen (1986) we need to treat endpoints of the grid differently:
            if (j==0 || j==N-1){
                if (j==0){
                    Prob[i*N+j] = normalCDF((y_grid[0]-lambda*y_grid[i]+omega/2)/sigma_e);
                }
                else {
                    Prob[i*N+j] = 1-normalCDF((y_grid[N-1]-lambda*y_grid[i]-omega/2)/sigma_e);
                }
                
            // Points inside the grid:
            } else {
                Prob[i*N+j] = normalCDF((y_grid[j]-lambda*y_grid[i]+omega/2)/sigma_e)-normalCDF((y_grid[j]-lambda*y_grid[i]-omega/2)/sigma_e);
            }
        }
    }

    return 0;
}

