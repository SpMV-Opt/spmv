/*
 *   CSR(Compressed Sparse Row) format for sparse matrix.
 *
 */

#ifndef _CSR_H_
#define _CSR_H_

// M: rows of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, double *nz_vals, int *column_index, int *row_start,
         double *x, double *y) {
  double y0;
  int i, j;
  // loop over the rows of sparse matrix
  for (i = 0; i < M; ++i) {
    y0 = y[i];
    for (j = row_start[i]; j < row_start[i + 1]; ++j) {
      y0 += nz_vals[j] * x[column_index[j]];
    }
    y[i] = y0;
  }
}

#endif // _CSR_H_
