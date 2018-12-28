/*
 *   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
 *   and copies it to stdout.
 *
 *   Usage:  a.out [filename] r c > output
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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <set>
#include <sys/time.h>
#ifdef PAPI
#include <papi.h>
#endif

#include "spmv_opt.h"
#include "utils.h"

#define NUM_EVENTS 3

#ifdef PAPI
// handle papi test_fail
static void test_fail(const char *file, int line, const char *call, int retval);
#endif

int main(int argc, char *argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s [martix-market-filename] r c\n", argv[0]);
    exit(1);
  }

  // perf wall time
  struct timeval begin, end;
  double time_elapsed = 0.0;

#ifdef PAPI
  // papi performance numbers
  float real_time, proc_time, mflops;
  long long flpins;
  int retval;
  long long counters[NUM_EVENTS];
  int PAPI_events[] = {
      PAPI_TOT_CYC, // total cycles
      PAPI_L2_DCM,  // level 2 data cache misses
      PAPI_L2_DCA,  // level 2 data cache accesses
      // PAPI_FP_OPS
  };

  // init papi lib
  if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    fprintf(stdout, "PAPI_library_init Error!\n");
    exit(-1);
  }
#endif // PAPI

  // rows, columns, nz: number of non-zero elems
  // r: tiny block row size, c: tiny block column size
  int rows, columns, nz;
  int r = atoi(argv[2]);
  int c = atoi(argv[3]);

  // get matrix size from input matrix market file
  get_matrix_size(argv[1], rows, columns, nz);
  // ceil original matrix size
  int Rows = r * (int)ceil(rows / r);
  int Columns = c * (int)ceil(columns / c);
  int bm = (int)ceil(rows / r);
  int bn = (int)ceil(columns / c);

#ifdef DEBUG
  fprintf(stdout, "rows:%d columns:%d nz:%d Rows:%d Columns:%d\n", rows, columns, nz, Rows, Columns);
#endif // DEBUG

  // source vector x random generation
  double *x = (double *)malloc(Columns * sizeof(double));
  if (!x) {
    fprintf(stdout, "%s:%d, fail to malloc x!\n", __FILE__, __LINE__);
    exit(-1);
  }
  rand_gen(Columns, x);
  // destination vector y for output
  double *y = (double *)malloc(Rows * sizeof(double));
  if (!y) {
    fprintf(stdout, "%s:%d, fail to malloc y!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(y, y + Rows, 0.0);

  /// ========= Naive impl for result correctness verify =============
#ifdef RESULT_VERIFY
  // I: x-axis row index, J: y-axis column index, val: non-zero value
  int *I, *J;
  double *val;
  I = (int *)malloc(nz * sizeof(int));
  if (!I) {
    fprintf(stdout, "%s:%d, fail to malloc I!\n", __FILE__, __LINE__);
    exit(-1);
  }
  J = (int *)malloc(nz * sizeof(int));
  if (!J) {
    fprintf(stdout, "%s:%d, fail to malloc J!\n", __FILE__, __LINE__);
    exit(-1);
  }
  val = (double *)malloc(nz * sizeof(double));
  if (!val) {
    fprintf(stdout, "%s:%d, fail to malloc val!\n", __FILE__, __LINE__);
    exit(-1);
  }

  // get matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);
  // input sparse matrix A
  double *A = (double *)malloc(Rows * Columns * sizeof(double));
  if (!A) {
    fprintf(stdout, "%s:%d, fail to malloc A!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(A, A + Rows * Columns, 0.0);
  for (int i = 0; i < nz; ++i) {
    A[I[i] * Columns + J[i]] = val[i];
  }

  // naive implement to check the correctness of other optimizers
  double *result = (double *)malloc(Rows * sizeof(double));
  if (!result) {
    fprintf(stdout, "%s:%d, fail to malloc result!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(result, result + Rows, 0.0);
  // 1. naive implement
  // retval = PAPI_start_counters(PAPI_events, NUM_EVENTS);
  gettimeofday(&begin, NULL);
  naive(Rows, Columns, A, x, result);
  gettimeofday(&end, NULL);
  // retval = PAPI_read_counters(counters, NUM_EVENTS);
  time_elapsed = (double)(end.tv_sec - begin.tv_sec) +
                 (end.tv_usec - begin.tv_usec) / 1000000.0;
  printf("Naive impl wall time: %.3lf s \n", time_elapsed);
  // assert(retval == PAPI_OK);
  // printf("Naive impl perf: L2 cache miss %lld, L2 cache misses ratio: %.3lf
  // %, in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2] *
  // 100, counters[0]);
#endif // RESULT_VERIFY
  /// ==========================================================

  record_t *records = (record_t *)malloc(nz * sizeof(record_t));
  if (!records) {
    fprintf(stdout, "%s:%d, fail to malloc records!\n", __FILE__, __LINE__);
    exit(-1);
  }
  // parse matrix non-zero values from input matrix market file
  get_records(argv[1], records);

  /// ===================== BCSR impl ==========================
  std::fill(y, y + Rows, 0.0);
  // reorder the records with increase order by the row
  records_reorder_by_rows(nz, records);
#ifdef DEBUG
  fprintf(stdout, "After sort:\n");
  for (int i = 0; i < nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].r, records[i].c, records[i].val);
  }
#endif // DEBUG
  std::vector<std::pair<int, int>> block_idx;
  get_bcsr_block_idx(records, nz, r, c, block_idx);
  int block_nz = block_idx.size();
  // transform sparse matrix into bcsr format
  double *b_values = (double*)malloc(r * c * block_nz * sizeof(double));
  if (!b_values) {
    fprintf(stdout, "%s:%d, fail to malloc b_values!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(b_values, b_values + r * c * block_nz, 0.0);
  int *b_col_idx = (int*)malloc(block_nz * sizeof(int));
  if (!b_col_idx) {
    fprintf(stdout, "%s:%d, fail to malloc b_col_idx!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(b_col_idx, b_col_idx + block_nz, 0);
  int *b_row_start = (int*)malloc((bm + 1) * sizeof(int));
  if (!b_row_start) {
    fprintf(stdout, "%s:%d, fail to malloc b_row_start!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(b_row_start, b_row_start + (bm + 1), 0);
  cvt2bcsr(Rows, nz, r, c, records, block_nz, block_idx, b_values, b_col_idx, b_row_start);
#ifdef DEBUG
  fprintf(stdout, "=== Rows: %d, Columns: %d, nz: %d, r: %d, c: %d\n", Rows, Columns, nz, r, c);
#endif // DEBUG
  // 2. BCSR implement
  // retval = PAPI_start_counters(PAPI_events, NUM_EVENTS);
  gettimeofday(&begin, NULL);
  bcsr_2x2(bm, b_row_start, b_col_idx, b_values, x, y);
  gettimeofday(&end, NULL);
  free(b_row_start);
  free(b_col_idx);
  free(b_values);
  // retval = PAPI_read_counters(counters, NUM_EVENTS);
  time_elapsed = (double)(end.tv_sec - begin.tv_sec) +
                 (end.tv_usec - begin.tv_usec) / 1000000.0;
  printf("BCSR impl wall time: %.3lf s \n", time_elapsed);
  // assert(retval == PAPI_OK);
  // printf("BCSR impl perf: L2 cache miss %lld, L2 cache misses ratio: %.3lf %, in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2] * 100, counters[0]);
#ifdef RESULT_VERIFY
  // check BCSR correctness
  if (check(Rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < Rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
#endif // RESULT_VERIFY
  /// ===========================================================

  // memory release
#ifdef RESULT_VERIFY
  free(I);
  free(J);
  free(val);
  free(A);
  free(result);
#endif // RESULT_VERIFY
  free(records);
  free(x);
  free(y);
  return 0;
}

#ifdef PAPI
static void test_fail(const char *file, int line, const char *call,
                      int retval) {
  printf("%s\tFAILED\nLine # %d\n", file, line);
  if (retval == PAPI_ESYS) {
    char buf[128];
    memset(buf, '\0', sizeof(buf));
    sprintf(buf, "System error in %s:", call);
    perror(buf);
  } else if (retval > 0) {
    printf("Error calculating: %s\n", call);
  } else {
    printf("Error in %s: %s\n", call, PAPI_strerror(retval));
  }
  printf("\n");
  exit(1);
}
#endif
