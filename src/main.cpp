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

#include "utils.h"
#include "spmv_opt.h"

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
#endif

  I = (int *)malloc(nz * sizeof(int));
  J = (int *)malloc(nz * sizeof(int));
  val = (double *)malloc(nz * sizeof(double));

  // parse matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);

  // source vector x random generation
  double *x = (double*) malloc(columns*sizeof(double));
  rand_gen(columns, x);
  // destination vector y for output
  double *y = (double*) malloc(rows*sizeof(double));
  std::fill(y, y + rows, 0.0);

#ifdef DEBUG
  for (int i = 0; i < nz; i++)
    fprintf(stdout, "%d %d %lf\n", I[i], J[i], val[i]);
#endif

#ifdef NAIVE
  // naive implement
  double *A = (double*)malloc(rows * columns * sizeof(double));
  std::fill(A, A + rows * columns, 0.0);
  for(int i = 0; i < nz; ++i) {
    A[I[i] * columns + J[i]] = val[i];
  }
  naive(rows, columns, A, x, y);
  // FIXME: check output correctness
#endif

  return 0;
}
