#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int read_matrix_market_to_csr(const char *filename, int *nrows, int *ncols, int *nnz,
                              int **row_ptr, int **col_ind, double **values);

typedef struct {
    double *values;
    int *indptr;
    int *indices;
} CSRMatrix;

typedef struct {
    double *values;
    int *indptr;
    int *indices;
} CSCMatrix;

CSCMatrix initialize_csc_matrix() {
    CSCMatrix matrix;
    matrix.values = NULL;
    matrix.indices = NULL;
    matrix.indptr = NULL;
    return matrix;
}
CSRMatrix initialize_csr_matrix() {
    CSRMatrix matrix;
    matrix.values = NULL;
    matrix.indices = NULL;
    matrix.indptr = NULL;
    return matrix;
}

static const CSCMatrix ERROR_CSC_MATRIX = {NULL, NULL, NULL};
static const CSRMatrix ERROR_CSR_MATRIX = {NULL, NULL, NULL};

CSRMatrix sparse_matrix_multiplication(CSRMatrix matrix, int *vector, int c, int n) {
    int s;
    double el;
    CSCMatrix csc = initialize_csc_matrix();
    int val, start_idx, end_idx;
    int max_nnz = matrix.indptr[n];
    csc.values = (double *)malloc(max_nnz * sizeof(double));
    csc.indices = (int *)malloc(max_nnz * sizeof(int));
    csc.indptr = (int *)malloc((c + 1) * sizeof(int));
    int nnz = 0;

    if (csc.values == NULL || csc.indices == NULL || csc.indptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(csc.values);
        free(csc.indices);
        free(csc.indptr);
        return ERROR_CSR_MATRIX;
    }

    for (int col = 0; col < c; col++) {
        csc.indptr[col] = nnz;  // Update column pointer

        for (int i = 0; i < n; i++) {
            val = 0;
            start_idx = matrix.indptr[i];
            end_idx = matrix.indptr[i + 1];

            for (int j = start_idx; j < end_idx; j++) {
                s = matrix.indices[j];
                el = matrix.values[j];

                if (vector[s] == col) {
                    val += el;
                }
            }

            if (val != 0) {
                csc.values[nnz] = val;
                csc.indices[nnz] = i;
                nnz++;
            }

        }
    }

    csc.indptr[c] = nnz;  // Update the last column pointer

    // Reallocate memory to the correct size
    csc.values = realloc(csc.values, nnz * sizeof(double));
    csc.indices = realloc(csc.indices, nnz * sizeof(int));

	//2os ypologismos

    s=0;
    el=0.0;
    CSRMatrix csr=initialize_csr_matrix();
    int nnz_res=0;
    max_nnz = nnz;  //to megisto nnz pou mporei na exei einai osa nnz eixe o csr prin.
    csr.values=(double *)malloc(max_nnz * sizeof(double));
    csr.indices=(int *)malloc(max_nnz * sizeof(int));
    csr.indptr=(int *)malloc((c + 1) * sizeof(int));
    if (csr.values == NULL || csr.indices == NULL || csr.indptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");

        free(csr.values);
        free(csr.indices);
        free(csr.indptr);
        return ERROR_CSR_MATRIX;
    }

    for (int row = 0; row < c; row++){
        csr.indptr[row] = nnz_res;
        for (int i = 0; i < c; i++){
            val=0;
            start_idx= csc.indptr[i];
            end_idx= csc.indptr[i + 1];
            for (int j = start_idx; j <end_idx; j++) {
                s = csc.indices[j];
                el=csc.values[j];
                if (vector[s] == row){
                    val+=el;
                }

            }
            if (val!=0){
                csr.values[nnz_res] = val;
                csr.indices[nnz_res] = i;
                nnz_res++;

            }


        }

    }
    csr.indptr[c]=nnz_res;

    csr.values = realloc(csr.values, nnz_res * sizeof(double));
    csr.indices = realloc(csr.indices, nnz_res * sizeof(int));
    free(csc.values);
    free(csc.indices);
    free(csc.indptr);
    return csr;
}

void print_csr_matrix(CSRMatrix matrix, int c) {
    printf("Result in CSR form:\n");


    printf("Values: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%f ", matrix.values[i]);
    }
    printf("\n");


    printf("Indices: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%d ", matrix.indices[i]);
    }
    printf("\n");


    printf("Indptr: ");
    for (int i = 0; i <= c; i++) {
        printf("%d ", matrix.indptr[i]);
    }
    printf("\n");
}
void print_csc_matrix(CSCMatrix matrix, int c) {
    printf("Result in CSC form:\n");


    printf("Values: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%f ", matrix.values[i]);
    }
    printf("\n");


    printf("Indices: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%d ", matrix.indices[i]);
    }
    printf("\n");


    printf("Indptr: ");
    for (int i = 0; i <= c; i++) {
        printf("%d ", matrix.indptr[i]);
    }
    printf("\n");
}

int main() {


    clock_t start_time = clock();
    int nrows, ncols, nnz;
    // CSR arrays
    int *row_ptr, *col_ind;
    double *values;

    //PARALLEL
    if (read_matrix_market_to_csr("af23560.mtx", &nrows, &ncols, &nnz, &row_ptr, &col_ind, &values) != 0) {
        return 1;  // Error reading matrix
    }

    // Now you can use the sparse matrix in CSR format for your functions
    CSRMatrix csr_matrix,csr_result;
    csr_matrix.values = values;
    csr_matrix.indices = col_ind;
    csr_matrix.indptr = row_ptr;

    // Create a vector with random values between 0 and ncols - 1
    int FIXED_SEED=5;
    srand(FIXED_SEED); // Seed for random number generation
    int *vector = (int *)malloc(nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
        vector[i] = rand()%ncols;
    }

    // Call the function from matrix_ops.c
    csr_result = sparse_matrix_multiplication(csr_matrix, vector, ncols, nrows);
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f seconds\n", elapsed_time);

    // Print the result
    //print_csc_matrix(csc_result, ncols);

    // Free allocated memory for CSR format
    free(row_ptr);
    free(col_ind);
    free(values);

    // Free the vector
    free(vector);

    // Free the allocated memory for the CSC matrix
    free(csr_result.values);
    free(csr_result.indices);
    free(csr_result.indptr);

    return 0;

}
/*
int main() {
    int num_rows = 40;  // Adjusted for a larger matrix
    int num_cols = 4;

    // Assuming dense_array is given in CSR form
    CSRMatrix csr_matrix,csr_result;
    int nnz = 18;  // Adjusted for a larger matrix
    csr_matrix.values = (double *)malloc(nnz * sizeof(double));
    csr_matrix.indices = (int *)malloc(nnz * sizeof(int));
    csr_matrix.indptr = (int *)malloc((num_rows + 1) * sizeof(int));


    double values[] = {1, 1, 2, 1, 1, 3, 1, 5, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int indices[] = {1, 2, 3, 0, 3, 2, 3, 0, 2, 1, 2, 3, 0, 1, 2, 3, 0, 1};
    int indptr[] = {0, 3, 5, 7, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9};

    for (int i = 0; i < nnz; i++) {
        csr_matrix.values[i] = values[i];
        csr_matrix.indices[i] = indices[i];
    }

    for (int i = 0; i < num_rows + 1; i++) {
        csr_matrix.indptr[i] = indptr[i];
    }

    // Create a vector
    int vector[] = {1, 1, 0, 0};

    // Call sparse_matrix_multiplication function directly with CSR matrix
    csr_result = sparse_matrix_multiplication(csr_matrix, vector, num_cols, num_rows);

    // Print the result
    print_csr_matrix(csr_result, num_cols);

    // Free the allocated memory for the CSC matrix
    free(csr_result.values);
    free(csr_result.indices);
    free(csr_result.indptr);

    // Free the allocated memory
    free(csr_matrix.values);
    free(csr_matrix.indices);
    free(csr_matrix.indptr);

    return 0;
}
*/
