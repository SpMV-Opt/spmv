/*
 *   Matrix Market I/O example program
 *
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.  This porgram does nothing useful, but
 *   illustrates common usage of the Matrix Matrix I/O routines.
 *   (See http://math.nist.gov/MatrixMarket for details.)
 *
 *   Usage:  a.out [filename] > output
 *
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

#include "utils.h"

#include <cstdlib>
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }

  int M, N, nz;
  int *I, *J;
  double *val;
  get_matrix_size(argv[1], M, N, nz);
  printf("%d %d %d\n", M, N, nz);

  I = (int *)malloc(nz * sizeof(int));
  J = (int *)malloc(nz * sizeof(int));
  val = (double *)malloc(nz * sizeof(double));

  get_matrix(argv[1], I, J, val);

  int i;
  for (i = 0; i < nz; i++)
    fprintf(stdout, "%d %d %lf\n", I[i], J[i], val[i]);

  return 0;
}
