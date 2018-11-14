/*
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.
 *
 *   Usage:  a.out [filename] > output
 *
 *   NOTES:
 *
 *   1) Matrix Market files are always 1-based, i.e. the index of the first
 *      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
 *      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
 *      to files.
 *
 *   2) ANSI C requires one to use the "l" format modifier when reading
 *      double precision floating point numbers in scanf() and
 *      its variants.  For example, use "%lf", "%lg", or "%le"
 *      when reading doubles, otherwise errors will occur.
 */

#include <cstdio>
#include <cstdlib>

#include "spmv_opt.h"
#include "utils.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }

  // rows, columns, nz: number of non-zero elems
  int rows, columns, nz;
  // I: x-axis row index, J: y-axis column index, val: non-zero value
  int *I, *J;
  float *val;

  // parse matrix size from input matrix market file
  get_matrix_size(argv[1], rows, columns, nz);
#ifdef DEBUG
  fprintf(stdout, "%d %d %d\n", rows, columns, nz);
#endif // DEBUG

  I = (int *)malloc(nz * sizeof(int));
  J = (int *)malloc(nz * sizeof(int));
  val = (float *)malloc(nz * sizeof(float));

  // parse matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);

  record_t *records = (record_t*)malloc(nz * sizeof(record_t));
  if(!records)
    fprintf(stdout, "ERROR! Fail to malloc records!\n");
  // parse matrix non-zero values from input matrix market file
  get_records(argv[1], records);
  // reorder the records with increase order by the row
  records_reorder_by_rows(nz, records);
  fprintf(stdout, "After sort:\n");
  for(int i = 0; i< nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].i, records[i].j, records[i].nz_val);
  }

  // source vector x random generation
  float *x = (float *)malloc(columns * sizeof(float));
  rand_gen(columns, x);
  // destination vector y for output
  float *y = (float *)malloc(rows * sizeof(float));
  std::fill(y, y + rows, 0.0);

  // input sparse matrix A
  float *A = (float *)malloc(rows * columns * sizeof(float));
  std::fill(A, A + rows * columns, 0.0);
  for (int i = 0; i < nz; ++i) {
    A[I[i] * columns + J[i]] = val[i];
  }

  // naive implement to check the correctness of other optimizers
  float *result = (float *)malloc(rows * sizeof(float));
  std::fill(result, result + rows, 0.0);
  naive(rows, columns, A, x, result);

#ifdef NAIVE
  // 1. naive implement
  naive(rows, columns, A, x, y);
  // FIXME: check output correctness
#endif // NAIVE

#ifdef CSR
  // FIXME: we assume each row must contain at least one element
  // transform sparse matrix into csr format
  float *nz_vals = (float *)malloc(nz * sizeof(float));
  int *column_index = (int *)malloc(nz * sizeof(int));
  int *row_start = (int *)malloc((rows + 1) * sizeof(int));
  //cvt2csr(rows, columns, nz, A, nz_vals, column_index, row_start);
  cvt2csr(rows, columns, nz, records, nz_vals, column_index, row_start);

#ifdef DEBUG
  fprintf(stdout, "rows: %d, columns: %d\n", rows, columns);
  // row_start[rows] = rows;
  int i;
  fprintf(stdout, "non-zero values:\n");
  for (i = 0; i < nz; ++i) {
    fprintf(stdout, "%lf ", nz_vals[i]);
  }
  fprintf(stdout, "\ncolumn index:\n");
  for (i = 0; i < nz; ++i) {
    fprintf(stdout, "%d ", column_index[i]);
  }
  fprintf(stdout, "\nrow start:\n");
  for (i = 0; i < rows + 1; ++i) {
    fprintf(stdout, "%d ", row_start[i]);
  }
  fprintf(stdout, "\n");
#endif // DEBUG
  // 2. CSR implement
  csr(rows, nz_vals, column_index, row_start, x, y);
  // check CSR correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
  free(nz_vals);
  free(column_index);
  free(row_start);
#endif // CSR

#ifdef CSC
  // transform sparse matrix into csc format
  float *nz_vals = (float *)malloc(nz * sizeof(float));
  int *row_index = (int *)malloc(nz * sizeof(int));
  int *column_start = (int *)malloc((columns + 1) * sizeof(int));
  cvt2csc(rows, columns, nz, A, nz_vals, row_index, column_start);
  // 3. CSC implement
  csc(columns, nz_vals, row_index, column_start, x, y);
  // check CSC correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
  free(nz_vals);
  free(row_index);
  free(column_start);

#endif // CSC

  // memory release
  free(I);
  free(J);
  free(val);
  free(A);
  free(records);
  free(x);
  free(y);
  free(result);
  return 0;
}
