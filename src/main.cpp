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
#include <cstring>
#include <cassert>
#include <papi.h>
#include <sys/time.h>

#include "spmv_opt.h"
#include "utils.h"

#define NUM_EVENTS 3

// handle papi test_fail 
static void test_fail(const char *file, int line, const char *call, int retval);

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
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
    //PAPI_FP_OPS
  };

  // init papi lib
  if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    fprintf(stdout, "PAPI_library_init Error!\n");
    exit(-1);
  }
#endif // PAPI

  // rows, columns, nz: number of non-zero elems
  int rows, columns, nz;

  // parse matrix size from input matrix market file
  get_matrix_size(argv[1], rows, columns, nz);

#ifdef DEBUG
  fprintf(stdout, "%d %d %d\n", rows, columns, nz);
#endif // DEBUG

  // source vector x random generation
  double *x = (double *)malloc(columns * sizeof(double));
  if (!x) {
    fprintf(stdout, "%s:%d, fail to malloc x!\n", __FILE__, __LINE__);
    exit(-1);
  }
  rand_gen(columns, x);
  // destination vector y for output
  double *y = (double *)malloc(rows * sizeof(double));
  if (!y) {
    fprintf(stdout, "%s:%d, fail to malloc y!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(y, y + rows, 0.0);

  /// =============== Naive impl for result verify ===================
#ifdef RESULT_VERIFY
  // I: x-axis row index, J: y-axis column index, val: non-zero value
  int *I, *J;
  double *val;
  I = (int *)malloc(nz * sizeof(int));
  if(!I) {
    fprintf(stdout, "%s:%d, fail to malloc I!\n", __FILE__, __LINE__);
    exit(-1);
  }
  J = (int *)malloc(nz * sizeof(int));
  if(!J){
    fprintf(stdout, "%s:%d, fail to malloc J!\n", __FILE__, __LINE__);
    exit(-1);
  }
  val = (double *)malloc(nz * sizeof(double));
  if(!val){
    fprintf(stdout, "%s:%d, fail to malloc val!\n", __FILE__, __LINE__);
    exit(-1);
  }

  // parse matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);
  // input sparse matrix A
  double *A = (double *)malloc(rows * columns * sizeof(double));
  if(!A) {
    fprintf(stdout, "%s:%d, fail to malloc A!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(A, A + rows * columns, 0.0);
  for (int i = 0; i < nz; ++i) {
    A[I[i] * columns + J[i]] = val[i];
  }

  // naive implement to check the correctness of other optimizers
  double *result = (double *)malloc(rows * sizeof(double));
  if(!result) {
    fprintf(stdout, "%s:%d, fail to malloc result!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(result, result + rows, 0.0);
  // 1. naive implement
  //retval = PAPI_start_counters(PAPI_events, NUM_EVENTS);
  gettimeofday(&begin, NULL);
  naive(rows, columns, A, x, result);
  gettimeofday(&end, NULL);
  //retval = PAPI_read_counters(counters, NUM_EVENTS);
  time_elapsed = (double)(end.tv_sec - begin.tv_sec)+ (end.tv_usec - begin.tv_usec)/1000000.0; 
  printf("Naive impl wall time: %.3lf ms \n", time_elapsed);
  //assert(retval == PAPI_OK);
  //printf("Naive impl perf: L2 cache miss %lld, L2 cache misses ratio: %.3lf %, in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2] * 100, counters[0]); 
#endif // RESULT_VERIFY
  /// ==========================================================

  record_t *records = (record_t *)malloc(nz * sizeof(record_t));
  if (!records) {
    fprintf(stdout, "%s:%d, fail to malloc records!\n", __FILE__, __LINE__);
    exit(-1);
  }
  // parse matrix non-zero values from input matrix market file
  get_records(argv[1], records);

  /// ================ CSR / CSR OMP impl ======================
  std::fill(y, y + rows, 0.0);
#if (defined CSR) || (defined CSR_OMP)
  // reorder the records with increase order by the row
  records_reorder_by_rows(nz, records);
#ifdef DEBUG
  fprintf(stdout, "After sort:\n");
  for (int i = 0; i < nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].r, records[i].c, records[i].val);
  }
#endif // DEBUG
  // transform sparse matrix into csr format
  double *nz_vals = (double *)malloc(nz * sizeof(double));
  if (!nz_vals) {
    fprintf(stdout, "%s:%d, fail to malloc nz_vals!\n", __FILE__, __LINE__);
    exit(-1);
  }
  int *column_index = (int *)malloc(nz * sizeof(int));
  if (!column_index) {
    fprintf(stdout, "%s:%d, fail to malloc column_index!\n", __FILE__,
            __LINE__);
    exit(-1);
  }
  int *row_start = (int *)malloc((rows + 1) * sizeof(int));
  if (!row_start) {
    fprintf(stdout, "%s:%d, fail to malloc row_start!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(row_start, row_start + (rows + 1), 0);
  cvt2csr(rows, nz, records, nz_vals, column_index, row_start);

#ifdef DEBUG
  fprintf(stdout, "rows: %d, columns: %d\n", rows, columns);
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
  //retval = PAPI_start_counters(PAPI_events, NUM_EVENTS);
  gettimeofday(&begin, NULL);
  csr(rows, nz_vals, column_index, row_start, x, y);
  gettimeofday(&end, NULL);
  //retval = PAPI_read_counters(counters, NUM_EVENTS);
  time_elapsed = (double)(end.tv_sec - begin.tv_sec)+ (end.tv_usec - begin.tv_usec)/1000000.0; 
  printf("CSR / CSR OMP impl wall time: %.3lf ms \n", time_elapsed);
  //assert(retval == PAPI_OK);
  //printf("CSR / CSR OMP impl perf: L2 cache miss %lld, L2 cache misses ratio: %.3lf %, in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2] * 100, counters[0]); 
#ifdef RESULT_VERIFY
  // check CSR correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
#endif // RESULT_VERIFY
  free(nz_vals);
  free(column_index);
  free(row_start);
#endif // CSR || CSR OMP
  /// ===========================================================

  /// ================ CSC / CSC OMP impl =======================
  std::fill(y, y + rows, 0.0);
#if (defined CSC) || (defined CSC_OMP)
  // reorder the records with increase order by the row
  records_reorder_by_columns(nz, records);
#ifdef DEBUG
  fprintf(stdout, "After sort:\n");
  for (int i = 0; i < nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].r, records[i].c, records[i].val);
  }
#endif // DEBUG
  // transform sparse matrix into csc format
  double *nz_vals = (double *)malloc(nz * sizeof(double));
  if (!nz_vals) {
    fprintf(stdout, "%s:%d, fail to malloc nz_vals!\n", __FILE__, __LINE__);
    exit(-1);
  }
  int *row_index = (int *)malloc(nz * sizeof(int));
  if (!row_index) {
    fprintf(stdout, "%s:%d, fail to malloc row_index!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(row_index, row_index + nz, 0);
  int *column_start = (int *)malloc((columns + 1) * sizeof(int));
  if (!column_start) {
    fprintf(stdout, "%s:%d, fail to malloc column_start!\n", __FILE__,
            __LINE__);
    exit(-1);
  }
  std::fill(column_start, column_start + (columns + 1), 0);
  cvt2csc(columns, nz, records, nz_vals, row_index, column_start);
#ifdef DEBUG
  fprintf(stdout, "rows: %d, columns: %d\n", rows, columns);
  int i;
  fprintf(stdout, "non-zero values:\n");
  for (i = 0; i < nz; ++i) {
    fprintf(stdout, "%lf ", nz_vals[i]);
  }
  fprintf(stdout, "\nrow index:\n");
  for (i = 0; i < nz; ++i) {
    fprintf(stdout, "%d ", row_index[i]);
  }
  fprintf(stdout, "\ncolumn start:\n");
  for (i = 0; i < columns + 1; ++i) {
    fprintf(stdout, "%d ", column_start[i]);
  }
  fprintf(stdout, "\n");
#endif // DEBUG

  // 3. CSC implement
  //retval = PAPI_start_counters(PAPI_events, NUM_EVENTS);
  gettimeofday(&begin, NULL);
  csc(columns, nz_vals, row_index, column_start, x, y);
  gettimeofday(&end, NULL);
  //retval = PAPI_read_counters(counters, NUM_EVENTS);
  //assert(retval == PAPI_OK);
  time_elapsed = (double)(end.tv_sec - begin.tv_sec)+ (end.tv_usec - begin.tv_usec)/1000000.0; 
  printf("CSC / CSC OMP impl wall time: %.3lf ms \n", time_elapsed);
  //printf("CSC / CSC OMP impl perf: L2 cache miss %lld, L2 cache misses ratio: %.3lf %, in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2] * 100, counters[0]); 
#ifdef RESULT_VERIFY
  // check CSC correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
#endif // RESULT_VERIFY
  free(nz_vals);
  free(row_index);
  free(column_start);
#endif // CSC || CSC OMP
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

static void test_fail(const char *file, int line, const char *call, int retval) {
  printf("%s\tFAILED\nLine # %d\n", file, line);
  if ( retval == PAPI_ESYS ) {
    char buf[128];
    memset( buf, '\0', sizeof(buf) );
    sprintf(buf, "System error in %s:", call );
    perror(buf);
  }
  else if ( retval > 0 ) {
    printf("Error calculating: %s\n", call );
  }
  else {
    printf("Error in %s: %s\n", call, PAPI_strerror(retval) );
  }
  printf("\n");
  exit(1);
}
