#ifndef JACOBI_H
#define JACOBI_H

#define single float
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>
// HPC Libraries
#include "/Users/cyae/Downloads/lapack-3.12.1/build/include/cblas.h"
#include "/Users/cyae/Downloads/lapack-3.12.1/build/include/lapacke.h"

#ifdef __cplusplus
extern "C" {
#endif
    void dlarge_(lapack_int *n, double *a, lapack_int *lda, lapack_int *iseed, double *work, lapack_int *info);
    void dlaror_ (char *side, char *init, lapack_int *m, lapack_int *n, double *a, lapack_int *lda, lapack_int *iseed, double *x, lapack_int *info);
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
        return Sigma;
    }

    vector<double> SSigma( m*n, 0.0);
    for (int i = 0; i < p; ++i) {
        SSigma[i + i * m] = sigma[i];
    }

    char side1 = 'L';
    char side2 = 'R';
    char init = 'N';
    vector<double> work2(2*m+n, 0.0);
    vector<double> work3(2*n+m, 0.0);
    dlaror_(&side1, &init, &m, &n, SSigma.data(), &m, iseed, work2.data(), &info);
    dlaror_(&side2, &init, &m, &n, SSigma.data(), &m, iseed, work3.data(), &info);

    return SSigma;
}

double cond2(std::vector<double> A, int m, int n) {
    int p = std::min(m, n);
    std::vector<double> s(p);
    std::vector<double> superb(std::max(1, p - 1));

    lapack_int info = LAPACKE_dgesvd(
            LAPACK_COL_MAJOR,
            'N', 'N',
            m, n,
            A.data(), m,
            s.data(),
            nullptr, 1,     // U not referenced when jobu='N'
            nullptr, 1,     // VT not referenced when jobvt='N'
            superb.data()   // must be allocated
    );

    if (info > 0) {
        std::cerr << "dgesvd did not converge, info = " << info << "\n";
        return std::numeric_limits<double>::infinity();
    }
    if (info < 0) {
        std::cerr << "dgesvd: illegal argument, info = " << info << "\n";
        return std::numeric_limits<double>::infinity();
    }
    if (s[p - 1] == 0.0) {return std::numeric_limits<double>::infinity();}
    return s[0] / s[p - 1];
}

single *right_singular_vector( const double* A, int m, int n);

#endif
