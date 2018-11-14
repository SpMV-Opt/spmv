/*
 *   Optimizers for sparse matrix multiply vector.
 *
 */

#ifndef _SPMV_OPT_H
#define _SPMV_OPT_H

// A: input sparse matrix, it has M x N float elements
// x: source vector, it has N x 1 float elements
// y: destination vector, it has M x 1 doule elements
void naive(const int &M, const int &N, float *A, float *x, float *y);

// M: rows of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 float elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, float *nz_vals, int *column_index, int *row_start,
         float *x, float *y);

// M: columns of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 float elements
// y: destination vector, it has M x 1 doule elements
void csc(const int &M, float *nz_vals, int *row_index, int *column_start,
         float *x, float *y);

#endif // _SPMV_OPT_H
