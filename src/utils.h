/*
 * Common utilities for matrix market file parser.
 *
 */
#ifndef _UTILS_H_
#define _UTILS_H_

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

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

#endif // _UTILS_H_
