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
  double *val;

  // parse matrix size from input matrix market file
  get_matrix_size(argv[1], rows, columns, nz);
#ifdef DEBUG
  fprintf(stdout, "%d %d %d\n", rows, columns, nz);
#endif // DEBUG

  I = (int *)malloc(nz * sizeof(int));
  J = (int *)malloc(nz * sizeof(int));
  val = (double *)malloc(nz * sizeof(double));

  // parse matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);

  // source vector x random generation
  double *x = (double *)malloc(columns * sizeof(double));
  rand_gen(columns, x);
  // destination vector y for output
  double *y = (double *)malloc(rows * sizeof(double));
  std::fill(y, y + rows, 0.0);

#ifdef DEBUG
  for (int i = 0; i < nz; i++)
    fprintf(stdout, "%d %d %lf\n", I[i], J[i], val[i]);
#endif // DEBUG

  // input sparse matrix A
  double *A = (double *)malloc(rows * columns * sizeof(double));
  std::fill(A, A + rows * columns, 0.0);
  for (int i = 0; i < nz; ++i) {
    A[I[i] * columns + J[i]] = val[i];
  }

#ifdef NAIVE
  // naive implement
  naive(rows, columns, A, x, y);
  // FIXME: check output correctness
#endif // NAIVE

#ifdef CSR
  // transform sparse matrix into csr format
  double *nz_vals = (double *)malloc(nz * sizeof(double));
  int *column_index = (int *)malloc(nz * sizeof(int));
  int *row_start = (int *)malloc((rows + 1) * sizeof(int));
  int count = 0;
  double element;
  for (int i = 0; i < rows; ++i) {
    row_start[i] = count;
    for (int j = 0; j < columns; ++j) {
      element = A[i * columns + j];
      if (element != 0.0) {
        column_index[count] = j;
        nz_vals[count] = element;
        ++count;
      }
    }
  }
  row_start[rows] = rows;
#ifdef DEBUG
  fprintf(stdout, "rows: %d, columns: %d\n", rows, columns);
  fprintf(stdout, "count: %d\n", count);
  //row_start[rows] = rows;
  int i;
  for(i = 0; i < nz; ++i) {
    fprintf(stdout, "%lf ", nz_vals[i]);
  }
  fprintf(stdout, "\n");
  for(i = 0; i < nz; ++i) {
    fprintf(stdout, "%d ", column_index[i]);
  }
  fprintf(stdout, "\n");
  for(i = 0; i < rows + 1; ++i) {
    fprintf(stdout, "%d ", row_start[i]);
  }
  fprintf(stdout, "\n");
#endif // DEBUG

  csr(rows, nz_vals, column_index, row_start, x, y);
#endif // CSR

  // memory release
  return 0;
}
