/*
 *   CSR(Compressed Sparse Row) format for sparse matrix.
 *
 */

// M: rows of sparse matrix
// column_block_size: column based block size
// nz_vals: non-zero elements of sparse matrix
// x: source vector, it has N x 1 float elements
// y: destination vector, it has M x 1 doule elements
void csr(const int &M, const int &column_block_size, float *nz_vals,
         int *column_index, int *row_start, float *x, float *y) {
  // loop over the column based block size
  for (int k = 0; k < column_block_size; ++k) {
    // perform csr spmv for each block
    float y0;
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
}
