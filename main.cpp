// CPP Own Libraries
#include <iostream>
#include <vector>

// HPC Libraries
#include "/Users/cyae/OpenBLAS/include/cblas.h"
#include "/Users/cyae/OpenBLAS/include/lapacke.h"
#include "jacobi.h"

using namespace std;

int main() {

    vector<double> a(10*10, 0.0);
    a= randsvd_work(10, 10, -1e10, 3);

    print_matrix(a.data(),10,10);


    return 0;
}
