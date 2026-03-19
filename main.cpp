// CPP Own Libraries
#include <iostream>
#include <vector>

// HPC Libraries
#include "/Users/cyae/OpenBLAS/include/cblas.h"
#include "/Users/cyae/OpenBLAS/include/lapacke.h"
#include "jacobi.h"

using namespace std;

int main() {


    std::vector<double> A = randsvd_work( 5, 5, -1e3, 3);

    print_matrix(A.data(), 5, 5);

    double cond = cond2( A, 5, 5);
    cout << cond << endl;

    return 0;

}
