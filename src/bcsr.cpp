/*
 *   Register Blocking CSR implement.
 *
 */

// bm: non-zero block rows of sparse matrix
void bcsr_1x1(int bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, double *x, double *y) {
  int i, j;
  double c0, d0;
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 1) {
      c0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * c0;
    }
    y[1 * i + 0] = d0;
  }
}

void bcsr_2x2(int bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, double *x, double *y) {
  int i, j;
  double d0, d1, c0, c1;
  /* loop over bm block rows */
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i];
    /* scalar replacement since reused */
    d1 = y[2 * i + 1];
    /* dense micro MVM */
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 2) {
      c0 = x[2 * b_col_idx[j] + 0]; /* scalar replacement since reused */
      c1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * c0;
      d1 += b_values[2] * c0;
      d0 += b_values[1] * c1;
      d1 += b_values[3] * c1;
    }
    y[2 * i] = d0;
    y[2 * i + 1] = d1;
  }
}

void bcsr_1x2(const int &bm, double *b_nz_vals, int *b_column_index,
              int *b_row_start, double *x, double *y) {
  int i, j;
  // loop over the block rows of sparse matrix
  for (i = 0; i < bm; ++i) {
    double t0 = 0;
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j) {
      t0 += b_nz_vals[j * 2 + 0] * x[b_column_index[j] + 0];
      t0 += b_nz_vals[j * 2 + 1] * x[b_column_index[j] + 1];
    }
    y[1 * i + 0] = t0;
  }
}

void bcsr_1x3(const int &bm, double *b_nz_vals, int *b_column_index,
              int *b_row_start, double *x, double *y) {
  int i, j;
  // loop over the block rows of sparse matrix
  for (i = 0; i < bm; ++i) {
    double t0 = 0;
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j) {
      t0 += b_nz_vals[j * 3 + 0] * x[b_column_index[j] + 0];
      t0 += b_nz_vals[j * 3 + 1] * x[b_column_index[j] + 1];
      t0 += b_nz_vals[j * 3 + 2] * x[b_column_index[j] + 2];
    }
    y[1 * i + 0] = t0;
  }
}
