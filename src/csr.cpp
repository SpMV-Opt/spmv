/*
 *   CSR(Compressed Sparse Row) format for sparse matrix.
 *
 */

// M: rows of sparse matrix
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 double elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, double *nz_vals, int *column_index, int *row_start,
         double *x, double *y) {
  double tmp;
  int i, j;
#ifdef CSR_OMP
#pragma omp parallel for default(shared) private(i, j, tmp)
#endif // CSR_OMP
  // loop over the rows of sparse matrix
  for (i = 0; i < M; ++i) {
    tmp = 0;
#ifdef CSR_OMP
    // hint compiler to ignore vector dependencies
    #pragma IVDEP
#endif // CSR_OMP
    for (j = row_start[i]; j < row_start[i + 1]; ++j) {
      tmp += nz_vals[j] * x[column_index[j]];
    }
    y[i] = tmp;
  }
}

