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

vector<double> randsvd_work(int m, int n, double kappa, int mode) {
    int posdef = 0;
    int p = min(m,n);
    vector<double> sigma(p, 0.0);

    if (kappa < 0){
        posdef = 1;
        kappa = -1.0*kappa;
    }

    switch (mode) {
        case 1:{
            fill(sigma.begin(), sigma.end(), 1/kappa);
            sigma[0] = 1.0;
            break;
        }
        case 2:{
            fill(sigma.begin(), sigma.end(), 1.0);
            sigma[p-1] = 1.0/kappa;
            break;
        }
        case 3:{
            double factor = pow(kappa, -1.0/(p-1));
//            fill(sigma.begin(), sigma.end(), factor);
            for (int i = 0; i < p; ++i) {sigma[i] = pow(factor, i);}
            break;
        }
        case 4:{
            fill(sigma.begin(), sigma.end(), 1.0);
            auto tmp = 1-1/kappa;
            for (int i = 0; i < p; i++){
                sigma[i] = sigma[i] - tmp*i/((p-1));
            }
            break;
        }
        default:{
            cout << "You need to specify the MODE.\n";
            break;
        }
    }

    vector<double> Sigma(p*p, 0.0);
    for (int i = 0; i < p; ++i) {
        Sigma[i + i * p] = sigma[i];
    }

    int iseed[4] = {1,2,3,5};
    vector<double> work(2*p, 0.0);
    int32_t info = 1;
    if (posdef){
        dlarge_(&p, Sigma.data(), &p, iseed, work.data(), &info);
    }
    return Sigma;
}

single *right_singular_vector( const double* A, int m, int n);

#endif
