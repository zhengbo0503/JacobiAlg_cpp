#!/bin/bash
set -e
cmake -B build
cmake --build build
echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
./build/one_sided_jacobi
