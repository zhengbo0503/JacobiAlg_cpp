#ifndef JACOBI_H
#define JACOBI_H

#define single float
#include <iostream>
#include <iomanip>
#include <vector>
// HPC Libraries
#include "/Users/cyae/Downloads/lapack-3.12.1/build/include/cblas.h"
#include "/Users/cyae/Downloads/lapack-3.12.1/build/include/lapacke.h"

#ifdef __cplusplus
extern "C" {
#endif
    void dlarge_(lapack_int *n, double *a, lapack_int *lda, lapack_int *iseed, double *work, lapack_int *info);
#ifdef __cplusplus
};
#endif

using namespace std;

template <typename T>
void print_matrix(const T* A, int m, int n) {
    cout << scientific << setprecision(2);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(12) << A[i + j * m] << " ";
        }
        cout << '\n';
    }
}

template<typename T>
T* randsvd_work(int m, int n, double kappa, int mode) {
    T* A = new T[m*n];
    int p, posdef;
    p = min(m,n);
//    vector<T> sigma(p,0);
    vector<T> ones(p, 1.0);


    if (kappa < 0){
        posdef = 1;
        kappa = -1.0*kappa;
    }

    switch (mode) {
        case 1:{
            vector<T> sigma(p, 1/kappa);
            sigma[1] = 1;
        }
        case 2:{
            vector<T> sigma(p, 1.0);
            sigma[p] = 1/kappa;
        }
        case 3:{
            double factor = pow(kappa, -1/(p-1));
            vector<T> sigma(p,factor);
            for (int i = 0; i < p; i++){ sigma[i] = pow(sigma[i], i); }
        }
        case 4:{
            vector<T> sigma(p,1.0);
            auto tmp = 1-1/kappa;
            for (int i = 0; i < p; i++){
                sigma[i] = sigma[i] - tmp*i/((p-1));
            }
        }
        default:{
            cout << "You need to specify the MODE.\n";
        }
    }



    return A;
}

single *right_singular_vector( const double* A, int m, int n);

#endif
