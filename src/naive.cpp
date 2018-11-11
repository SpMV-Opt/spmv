/*
 *   Naive implement for sparse matrix multiply vector.
 *
 */
#ifndef _NAIVE_H_
#define _NAIVE_H_

#include "spmv_opt.h"

// A: input sparse matrix, it has M x N double elements
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void naive(const int &M, const int &N, double *A, double *x, double *y) {
  for (int i = 0; i < M; ++i) {
    y[i] = 0;
    for (int j = 0; j < N; ++j) {
      y[i] += A[i * N + j] * x[j];
    }
  }
}

#endif // _NAIVE_H_
