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
#include <vector>

#define RANDOM_MAX 99
#define ESP 1e-6

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
