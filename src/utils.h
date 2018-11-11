/*
 * Common utilities for matrix market file parser.
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
void get_matrix(const char *file_name, int *I, int *J, double *val) {
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
    sscanf(line.c_str(), "%d %d %lf", &I[i], &J[i], &val[i]);
    I[i] -= 1;
    J[i] -= 1;
    ++i;
  }
  in.close();
}

// source vector x random generation
void rand_gen(const int &len, double *x) {
  std::mt19937 seed;
  std::uniform_real_distribution<> dis(-RAND_MAX, RAND_MAX);
  for (int i = 0; i < len; ++i) {
    x[i] = dis(seed);
    fprintf(stdout, "%lf ", x[i]);
  }
}

// check the correctness of optimizer output and naive result
bool check(const size_t &len, double *output, double *result) {
  bool pass = true;
  for(std::size_t i = 0; i < len; ++i) {
    if((output[i] - result[i]) > ESP) {
        pass = false;
        break;
    }
  }
  return pass;
}

#endif // _UTILS_H_
