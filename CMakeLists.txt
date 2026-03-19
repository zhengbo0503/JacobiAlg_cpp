cmake_minimum_required(VERSION 3.31)
project(one_sided_jacobi)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_COMPILER clang++)
set(OPENBLAS_INC "/Users/cyae/OpenBLAS/include")
set(OPENBLAS_LIB "/Users/cyae/OpenBLAS/lib/libopenblas.dylib")

add_executable(one_sided_jacobi main.cpp jacobi.h)

target_include_directories(one_sided_jacobi PRIVATE ${OPENBLAS_INC})
target_link_libraries(one_sided_jacobi PRIVATE ${OPENBLAS_LIB})
target_link_libraries(one_sided_jacobi PRIVATE openblas)
