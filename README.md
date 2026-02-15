# Sparse Matrix Multiplication (Sequential)
Overview

This C project implements efficient sparse matrix multiplication using CSR and CSC formats.
Intermediate matrices are computed indirectly via a vector that indicates the cluster assignment of each row. During multiplication, this vector is used to accumulate values directly into the intermediate sparse matrix, avoiding the creation of large sparse matrices, reducing memory usage, and improving cache efficiency.

This project was developed for the Parallel and Distributed Systems course at ECE AUTH (2023–2024).

Input matrices are in Matrix Market (.mtx) format.


Build & Run

Compile:

gcc sparse_mult_main.c sparse_utils.c mmio.c -o sparse_mult

Run:

./sparse_mult


Output:

Elapsed time: 149.137000 seconds

