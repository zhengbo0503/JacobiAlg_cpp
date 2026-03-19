#include "jacobi.h"
#include <iostream>
#include <iomanip>

#include "/Users/cyae/OpenBLAS/include/cblas.h"
#include "/Users/cyae/OpenBLAS/include/lapacke.h"

using namespace std;
#define single float






/*
 * lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
                           lapack_int m, lapack_int n, double* a,
                           lapack_int lda, double* s, double* u, lapack_int ldu,
                           double* vt, lapack_int ldvt, double* superb );
 */

single *right_singular_vector( const double* A, int m, int n) {

    int32_t size = m*n;
    single* As = new single[size];
    for (int i = 0; i < size; i++){As[i] = static_cast<single>(A[i]);}
    int32_t info_dgesvd = 0;
    single* s = new single[min(m,n)]();
    single* u = new single[m*m]();
    single* vt = new single[n*n]();
    single* superb = new single[min(m,n)-1]();
    info_dgesvd = LAPACKE_sgesvd(
            LAPACK_COL_MAJOR, 'N', 'A', m, n, As,
                                  'm', s, u, m, vt, n, superb);
    if (info_dgesvd < 0){
        cout << "ERROR : jacobi.cpp >> right_singular_vector >> LAPACKE_dgesvd\n";
        cout << "\tThe" << -1*info_dgesvd << "-th parameter had an illegal value.\n";
    } else if ( info_dgesvd > 0 ){
        cout << "ERROR : jacobi.cpp >> right_singular_vector >> LAPACKE_dgesvd\n";
        cout << "\tDBDSQR did not converge.\n";
    };
    return vt;
}

double *one_sided_jacobi(double *A) {
    // low precision svd

    // orthogonalization

    // precondition

    // jacobi
    return 0;
}
