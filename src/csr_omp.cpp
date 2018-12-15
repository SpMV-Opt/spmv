/*
 *   CSR(Compressed Sparse Row) format for sparse matrix.
 *
 */
#include <omp.h>

// M: rows of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, double *nz_vals, int *column_index, int *row_start,
         double *x, double *y) {
  double tmp;
  int i, j;
  #pragma omp parallel for default(shared) private(i, j, tmp)
  // loop over the rows of sparse matrix
  for (i = 0; i < M; ++i) {
    tmp = 0;
    for (j = row_start[i]; j < row_start[i + 1]; ++j) {
      tmp += nz_vals[j] * x[column_index[j]];
    }
    y[i] = tmp;
  }
}

