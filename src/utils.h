/*
 * Common utilities.
 *
 */
#ifndef _UTILS_H_
#define _UTILS_H_

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#define RANDOM_MAX 49
#define ESP 1e-6

typedef struct {
  int i, j;     // rows and columns
  float nz_val; // non-zero float value
} record_t;

bool cmp_key_i(const record_t &a, const record_t &b) { return a.i < b.i; }

bool cmp_key_j(const record_t &a, const record_t &b) { return a.j < b.j; }

// M: rows, N: columns, nz: number of non-zero elems
void get_matrix_size(const char *file_name, int &M, int &N, int &nz) {
  std::ifstream in;
  in.open(file_name);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", file_name);
  }
  std::string line;
  while (std::getline(in, line) && line != "") {
    // skip the comments
    if (line[0] == '%')
      continue;
    sscanf(line.c_str(), "%d %d %d", &M, &N, &nz);
    break;
  }
  in.close();
}

// nz: number of non-zero elems, I: x-axis, J: y-axis, val: value
void get_matrix(const char *file_name, int *I, int *J, float *val) {
  std::ifstream in;
  in.open(file_name);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", file_name);
  }
  int M, N, nz;
  std::string line;
  while (std::getline(in, line) && line != "") {
    // skip the comments
    if (line[0] == '%')
      continue;
    sscanf(line.c_str(), "%d %d %d", &M, &N, &nz);
    break;
  }

  int i = 0;
  while (std::getline(in, line) && line != "") {
    sscanf(line.c_str(), "%d %d %f", &I[i], &J[i], &val[i]);
    I[i] -= 1;
    J[i] -= 1;
    ++i;
  }
  in.close();
}

// nz: number of non-zero elems, I: x-axis, J: y-axis, val: value
void get_records(const char *file_name, record_t *records) {
  std::ifstream in;
  in.open(file_name);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", file_name);
  }
  int M, N, nz;
  std::string line;
  while (std::getline(in, line) && line != "") {
    // skip the comments
    if (line[0] == '%')
      continue;
    sscanf(line.c_str(), "%d %d %d", &M, &N, &nz);
    break;
  }

  int count = 0;
  record_t tmp;
  while (std::getline(in, line) && line != "") {
    sscanf(line.c_str(), "%d %d %f", &tmp.i, &tmp.j, &tmp.nz_val);
    tmp.i -= 1;
    tmp.j -= 1;
    records[count] = tmp;
    ++count;
  }
  in.close();
}

void records_reorder_by_rows(const int &nz, record_t *records) {
  std::sort(records, records + nz, cmp_key_i);
}

// source vector x random generation
void rand_gen(const int &len, float *x) {
  static thread_local std::mt19937 seed;
  std::uniform_real_distribution<> dis(-RAND_MAX, RAND_MAX);
  for (int i = 0; i < len; ++i) {
    x[i] = dis(seed);
    // fprintf(stdout, "%lf ", x[i]);
  }
}

// check the correctness of optimizer output and naive result
bool check(const size_t &len, float *output, float *result) {
  bool pass = true;
  for (std::size_t i = 0; i < len; ++i) {
    if ((output[i] - result[i]) > ESP) {
      pass = false;
      break;
    }
  }
  return pass;
}

void cvt2csr(const int &rows, const int &columns, const int &nz,
             record_t *records, float *nz_vals, int *column_index,
             int *row_start) {
  record_t *tmp = (record_t *)malloc(columns * sizeof(record_t));
  if (!tmp) {
    fprintf(stdout, "%s:%d fail to malloc tmp!\n", __FILE__, __LINE__);
    return;
  }
  int row_count = 1, row_offset = 0;
  int row_cur = records[0].i;
  row_start[0] = 0;
  tmp[0] = records[0];
  // FIXME: we assume nz >= 1
  for (int k = 1; k < nz; ++k) {
    if (row_cur == records[k].i) {
      tmp[row_count] = records[k];
      row_cur = records[k].i;
      ++row_count;
    } else {
      // get row_start
      ++row_offset;
      row_start[row_offset] = row_start[row_offset - 1] + row_count;
      // sort column index by j
      std::sort(tmp, tmp + row_count, cmp_key_j);
      // get column_index and nz_vals
      for (int m = 0; m < row_count; ++m) {
        column_index[k - row_count + m] = tmp[m].j;
        nz_vals[k - row_count + m] = tmp[m].nz_val;
      } // end for
      row_count = 1;
    }   // end else
  }     // end for
  row_start[rows] = nz;
  free(tmp);
}

#if 0
// convert the input matrix into csr format
void cvt2csr(const int &rows, const int &columns, const int &nz, float *A,
             float *nz_vals, int *column_index, int *row_start) {
  int count = 0;
  float element;
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
  row_start[rows] = nz;
}
#endif

// convert the input matrix into csc format
void cvt2csc(const int &rows, const int &columns, const int &nz, float *A,
             float *nz_vals, int *row_index, int *column_start) {
  int count = 0;
  float element;
  for (int j = 0; j < columns; ++j) {
    column_start[j] = count;
    for (int i = 0; i < rows; ++i) {
      element = A[i * columns + j];
      if (element != 0.0) {
        row_index[count] = i;
        nz_vals[count] = element;
        ++count;
      }
    }
  }
  column_start[columns] = nz;
}

#endif // _UTILS_H_
