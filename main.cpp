// CPP Own Libraries
#include <iostream>
#include <vector>

// HPC Libraries
#include "/Users/cyae/OpenBLAS/include/cblas.h"
#include "/Users/cyae/OpenBLAS/include/lapacke.h"
#include "jacobi.h"

using namespace std;

int main() {
    vector<double> ones(20, 1.0);

    int n = 5;
    int lda = n;
    int info = 1;
    int iseed[4] = {1,2,3,5};
    double* A = new double[lda*n];
    double* work = new double[2 * n];
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            A[i + j * lda] = (i == j) ? (i + 1.0) : 0.0;
        }
    }
    print_matrix(A, n, lda);

    dlarge_(&n, A, &lda, iseed, work, &info);

    print_matrix(A, n, lda);
    cout << info << endl;

    return 0;
}
