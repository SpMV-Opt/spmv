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

#include "spmv_opt.h"
#include "utils.h"
#include "papi.h"

//handle papi test_fail 
static void test_fail(const char *file, int line, const char *call, int retval);

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }
  
  //papi performance numbers 
  float real_time, proc_time, mflops; 
  long long flpins;
  int retval;

  long long counters[4];
  int PAPI_events[] = {
      PAPI_TOT_CYC, //total cycles
      PAPI_L2_DCM,  //level 2 data cache misses
      PAPI_L2_DCA,   //level 2 data cache accesses
    //  PAPI_FP_OPS
  };
  // rows, columns, nz: number of non-zero elems
  int rows, columns, nz;
#if 0
  // I: x-axis row index, J: y-axis column index, val: non-zero value
  int *I, *J;
#endif
  double *val;
  
  //init papi lib
  PAPI_library_init(PAPI_VER_CURRENT);

  // parse matrix size from input matrix market file
  get_matrix_size(argv[1], rows, columns, nz);
#ifdef DEBUG
  fprintf(stdout, "%d %d %d\n", rows, columns, nz);
#endif // DEBUG

#if 0
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
  val = (float *)malloc(nz * sizeof(float));
  if(!val){
    fprintf(stdout, "%s:%d, fail to malloc val!\n", __FILE__, __LINE__);
    exit(-1);
  }

  // parse matrix non-zero values from input matrix market file
  get_matrix(argv[1], I, J, val);
#endif

  record_t *records = (record_t *)malloc(nz * sizeof(record_t));
  if (!records) {
    fprintf(stdout, "%s:%d, fail to malloc records!\n", __FILE__, __LINE__);
    exit(-1);
  }
  // parse matrix non-zero values from input matrix market file
  get_records(argv[1], records);

  // source vector x random generation
  float *x = (float *)malloc(columns * sizeof(float));
  if (!x) {
    fprintf(stdout, "%s:%d, fail to malloc x!\n", __FILE__, __LINE__);
    exit(-1);
  }
  rand_gen(columns, x);
  // destination vector y for output
  float *y = (float *)malloc(rows * sizeof(float));
  if (!y) {
    fprintf(stdout, "%s:%d, fail to malloc y!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(y, y + rows, 0.0);

#if 0
  // input sparse matrix A
  float *A = (float *)malloc(rows * columns * sizeof(float));
  if(!A) {
    fprintf(stdout, "%s:%d, fail to malloc A!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(A, A + rows * columns, 0.0);
  double *y = (double *)malloc(rows * sizeof(double));
  std::fill_n(y, rows, 0.0);

#ifdef DEBUG
  //for (int i = 0; i < nz; i++)
  //  fprintf(stdout, "%d %d %lf\n", I[i], J[i], val[i]);
#endif // DEBUG

  for (int i = 0; i < nz; ++i) {
    A[I[i] * columns + J[i]] = val[i];
  }

  // naive implement to check the correctness of other optimizers
  float *result = (float *)malloc(rows * sizeof(float));
  if(!result) {
    fprintf(stdout, "%s:%d, fail to malloc result!\n", __FILE__, __LINE__);
    exit(-1);
  }
  std::fill(result, result + rows, 0.0);
  naive(rows, columns, A, x, result);
#endif

#ifdef NAIVE
#if 0
  // 1. naive implement
  naive(rows, columns, A, x, y);
#endif
 
  //Setup PAPI library and begin collecting data from the counters 
  //if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
  //  test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
 
  retval = PAPI_start_counters(PAPI_events, 3);
  // naive implement
  naive(rows, columns, A, x, y);
 
  retval = PAPI_read_counters(counters, 3);

  printf("%lld L2 cache misses (%.3lf%% misses) in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2], counters[0]); 
  //printf("FLOPS: %lld\n", counters[3]);
  // FIXME: check output correctness
  
  //Collect the data into the variables passed in */
  //if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
  //  test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
  
 // printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
 // real_time, proc_time, flpins, mflops);
 // printf("%s\tPASSED\n", __FILE__);
 // PAPI_shutdown();
#endif // NAIVE

#ifdef CSR
  // reorder the records with increase order by the row
  records_reorder_by_rows(nz, records);
  fprintf(stdout, "After sort:\n");
  for (int i = 0; i < nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].r, records[i].c, records[i].val);
  }
  // transform sparse matrix into csr format
  float *nz_vals = (float *)malloc(nz * sizeof(float));
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
  retval = PAPI_start_counters(PAPI_events, 3);
  // 2. CSR implement
  csr(rows, nz_vals, column_index, row_start, x, y);
  
  retval = PAPI_read_counters(counters, 3);
  printf("%lld L2 cache misses (%.3lf%% misses) in %lld cycles\n", counters[1], (double)counters[1]/(double)counters[2], counters[0]); 
#if 0
  // check CSR correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
#endif // 0
  free(nz_vals);
  free(column_index);
  free(row_start);
#endif // CSR

#ifdef CSC
  // reorder the records with increase order by the row
  records_reorder_by_columns(nz, records);
  fprintf(stdout, "After sort:\n");
  for (int i = 0; i < nz; ++i) {
    fprintf(stdout, "%d %d %f\n", records[i].r, records[i].c, records[i].val);
  }
  // transform sparse matrix into csc format
  float *nz_vals = (float *)malloc(nz * sizeof(float));
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
  csc(columns, nz_vals, row_index, column_start, x, y);
#if 0
  // check CSC correctness
  if (check(rows, y, result)) {
    fprintf(stdout, "PASS\n");
  } else {
    for (int i = 0; i < rows; ++i)
      fprintf(stdout, "y: %lf result: %lf\n", y[i], result[i]);
    fprintf(stdout, "FAILED\n");
  }
#endif
  free(nz_vals);
  free(row_index);
  free(column_start);

#endif // CSC

  // memory release
#if 0
  free(I);
  free(J);
  free(val);
  free(A);
  free(result);
#endif
  free(records);
  free(x);
  free(y);
  return 0;
}

static void test_fail(const char *file, int line, const char *call, int retval){
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
