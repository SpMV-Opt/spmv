/*
 *   Optimizers for sparse matrix multiply vector.
 *
 */

#ifndef _SPMV_OPT_H
#define _SPMV_OPT_H

typedef void (*BCSR_FP)(const int &bm, const int *b_row_start,
                        const int *b_col_idx, const double *b_values,
                        const double *x, double *y);
extern BCSR_FP bcsr_routines[12][12][1];

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

// M: columns of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void csc(const int &M, double *nz_vals, int *row_index, int *column_start,
         double *x, double *y);

#define DEF_BCSR_FUNC(SUFFIX)                                                  \
  void bcsr_##SUFFIX(const int &bm, const int *b_row_start,                    \
                     const int *b_col_idx, const double *b_values,             \
                     const double *x, double *y);

DEF_BCSR_FUNC(1x1)
DEF_BCSR_FUNC(1x2)
DEF_BCSR_FUNC(1x3)
DEF_BCSR_FUNC(1x4)
DEF_BCSR_FUNC(1x5)
DEF_BCSR_FUNC(1x6)
DEF_BCSR_FUNC(1x7)
DEF_BCSR_FUNC(1x8)
DEF_BCSR_FUNC(1x9)
DEF_BCSR_FUNC(1x10)
DEF_BCSR_FUNC(1x11)
DEF_BCSR_FUNC(1x12)
DEF_BCSR_FUNC(2x1)
DEF_BCSR_FUNC(2x2)
DEF_BCSR_FUNC(2x3)
DEF_BCSR_FUNC(2x4)
DEF_BCSR_FUNC(2x5)
DEF_BCSR_FUNC(2x6)
DEF_BCSR_FUNC(2x7)
DEF_BCSR_FUNC(2x8)
DEF_BCSR_FUNC(2x9)
DEF_BCSR_FUNC(2x10)
DEF_BCSR_FUNC(2x11)
DEF_BCSR_FUNC(2x12)
DEF_BCSR_FUNC(3x1)
DEF_BCSR_FUNC(3x2)
DEF_BCSR_FUNC(3x3)
DEF_BCSR_FUNC(3x4)
DEF_BCSR_FUNC(3x5)
DEF_BCSR_FUNC(3x6)
DEF_BCSR_FUNC(3x7)
DEF_BCSR_FUNC(3x8)
DEF_BCSR_FUNC(3x9)
DEF_BCSR_FUNC(3x10)
DEF_BCSR_FUNC(3x11)
DEF_BCSR_FUNC(3x12)
DEF_BCSR_FUNC(4x1)
DEF_BCSR_FUNC(4x2)
DEF_BCSR_FUNC(4x3)
DEF_BCSR_FUNC(4x4)
DEF_BCSR_FUNC(4x5)
DEF_BCSR_FUNC(4x6)
DEF_BCSR_FUNC(4x7)
DEF_BCSR_FUNC(4x8)
DEF_BCSR_FUNC(4x9)
DEF_BCSR_FUNC(4x10)
DEF_BCSR_FUNC(4x11)
DEF_BCSR_FUNC(4x12)
DEF_BCSR_FUNC(5x1)
DEF_BCSR_FUNC(5x2)
DEF_BCSR_FUNC(5x3)
DEF_BCSR_FUNC(5x4)
DEF_BCSR_FUNC(5x5)
DEF_BCSR_FUNC(5x6)
DEF_BCSR_FUNC(5x7)
DEF_BCSR_FUNC(5x8)
DEF_BCSR_FUNC(5x9)
DEF_BCSR_FUNC(5x10)
DEF_BCSR_FUNC(5x11)
DEF_BCSR_FUNC(5x12)
DEF_BCSR_FUNC(6x1)
DEF_BCSR_FUNC(6x2)
DEF_BCSR_FUNC(6x3)
DEF_BCSR_FUNC(6x4)
DEF_BCSR_FUNC(6x5)
DEF_BCSR_FUNC(6x6)
DEF_BCSR_FUNC(6x7)
DEF_BCSR_FUNC(6x8)
DEF_BCSR_FUNC(6x9)
DEF_BCSR_FUNC(6x10)
DEF_BCSR_FUNC(6x11)
DEF_BCSR_FUNC(6x12)
DEF_BCSR_FUNC(7x1)
DEF_BCSR_FUNC(7x2)
DEF_BCSR_FUNC(7x3)
DEF_BCSR_FUNC(7x4)
DEF_BCSR_FUNC(7x5)
DEF_BCSR_FUNC(7x6)
DEF_BCSR_FUNC(7x7)
DEF_BCSR_FUNC(7x8)
DEF_BCSR_FUNC(7x9)
DEF_BCSR_FUNC(7x10)
DEF_BCSR_FUNC(7x11)
DEF_BCSR_FUNC(7x12)
DEF_BCSR_FUNC(8x1)
DEF_BCSR_FUNC(8x2)
DEF_BCSR_FUNC(8x3)
DEF_BCSR_FUNC(8x4)
DEF_BCSR_FUNC(8x5)
DEF_BCSR_FUNC(8x6)
DEF_BCSR_FUNC(8x7)
DEF_BCSR_FUNC(8x8)
DEF_BCSR_FUNC(8x9)
DEF_BCSR_FUNC(8x10)
DEF_BCSR_FUNC(8x11)
DEF_BCSR_FUNC(8x12)
DEF_BCSR_FUNC(9x1)
DEF_BCSR_FUNC(9x2)
DEF_BCSR_FUNC(9x3)
DEF_BCSR_FUNC(9x4)
DEF_BCSR_FUNC(9x5)
DEF_BCSR_FUNC(9x6)
DEF_BCSR_FUNC(9x7)
DEF_BCSR_FUNC(9x8)
DEF_BCSR_FUNC(9x9)
DEF_BCSR_FUNC(9x10)
DEF_BCSR_FUNC(9x11)
DEF_BCSR_FUNC(9x12)
DEF_BCSR_FUNC(10x1)
DEF_BCSR_FUNC(10x2)
DEF_BCSR_FUNC(10x3)
DEF_BCSR_FUNC(10x4)
DEF_BCSR_FUNC(10x5)
DEF_BCSR_FUNC(10x6)
DEF_BCSR_FUNC(10x7)
DEF_BCSR_FUNC(10x8)
DEF_BCSR_FUNC(10x9)
DEF_BCSR_FUNC(10x10)
DEF_BCSR_FUNC(10x11)
DEF_BCSR_FUNC(10x12)
DEF_BCSR_FUNC(11x1)
DEF_BCSR_FUNC(11x2)
DEF_BCSR_FUNC(11x3)
DEF_BCSR_FUNC(11x4)
DEF_BCSR_FUNC(11x5)
DEF_BCSR_FUNC(11x6)
DEF_BCSR_FUNC(11x7)
DEF_BCSR_FUNC(11x8)
DEF_BCSR_FUNC(11x9)
DEF_BCSR_FUNC(11x10)
DEF_BCSR_FUNC(11x11)
DEF_BCSR_FUNC(11x12)
DEF_BCSR_FUNC(12x1)
DEF_BCSR_FUNC(12x2)
DEF_BCSR_FUNC(12x3)
DEF_BCSR_FUNC(12x4)
DEF_BCSR_FUNC(12x5)
DEF_BCSR_FUNC(12x6)
DEF_BCSR_FUNC(12x7)
DEF_BCSR_FUNC(12x8)
DEF_BCSR_FUNC(12x9)
DEF_BCSR_FUNC(12x10)
DEF_BCSR_FUNC(12x11)
DEF_BCSR_FUNC(12x12)

#endif // _SPMV_OPT_H
