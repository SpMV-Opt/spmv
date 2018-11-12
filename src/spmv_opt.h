/*
 *   Optimizers for sparse matrix multiply vector.
 *
 */

#ifndef _SPMV_OPT_H
#define _SPMV_OPT_H

// A: input sparse matrix, it has M x N double elements
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void naive(const int &M, const int &N, double *A, double *x, double *y);

// M: rows of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, double *nz_vals, int *column_index, int *row_start,
         double *x, double *y);


#endif // _SPMV_OPT_H
