/*
 * Common utilities.
 *
 */
#ifndef _UTILS_H_
#define _UTILS_H_

#include <algorithm>
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

#define RANDOM_MAX 2
#define ESP 1e-6

typedef struct {
  int r, c;   // rows and columns
  double val; // non-zero double value
} record_t;

typedef struct {
  int r, c; // rows and columns
} block_t;

bool cmp_record_key_row(const record_t &a, const record_t &b) {
  return a.r < b.r;
}

bool cmp_record_key_column(const record_t &a, const record_t &b) {
  return a.c < b.c;
}

bool cmp_record_row_col(const record_t &a, const record_t &b) {
  if (a.r != b.r)
    return a.r < b.r;
  return a.c < b.c;
}

bool cmp_block_key_row(const block_t &a, const block_t &b) { return a.r < b.r; }

bool cmp_block_key_column(const block_t &a, const block_t &b) {
  return a.c < b.c;
}

bool cmp_block_pair(const std::pair<int, int> &a,
                    const std::pair<int, int> &b) {
  if (a.first != b.first)
    return a.first < b.first;
  return a.second < b.second;
}

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

void get_bcsr_block_idx(const int &ROWS, const int &COLS, const char *file_name,
                        const int &nz, const int &_r, const int &_c,
                        std::vector<std::pair<int, int>> &block_idx) {
  int block_col = COLS / _c;
  int block_row = ROWS / _r;
  int bit_size = block_row * block_col;
  // std::bitset<bit_size> bit_map;
  bool *bit_map = (bool *)malloc(bit_size * sizeof(bool));
  std::fill(bit_map, bit_map + bit_size, false);
  if (!bit_map) {
    fprintf(stdout, "%s:%d, fail to malloc bit_map!\n", __FILE__, __LINE__);
    return;
  }
  std::ifstream in;
  in.open(file_name);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", file_name);
  }
  int _M, _N, _nz;
  std::string line;
  while (std::getline(in, line) && line != "") {
    // skip the comments
    if (line[0] == '%')
      continue;
    sscanf(line.c_str(), "%d %d %d", &_M, &_N, &_nz);
    break;
  }
  int _row, _col;
  double _val;
  for (int i = 0; i < nz; ++i) {
    std::getline(in, line);
    sscanf(line.c_str(), "%d %d %lf", &_row, &_col, &_val);
    _row -= 1;
    _col -= 1;
    bit_map[_row / _r * block_col + _col / _c] = true;
  }
  in.close();
  for (int i = 0; i < bit_size; ++i) {
    if (bit_map[i]) {
      block_idx.push_back({i / block_col, i % block_col});
    }
  }
  free(bit_map);
  std::sort(block_idx.begin(), block_idx.end(), cmp_block_pair);
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
    sscanf(line.c_str(), "%d %d %lf", &tmp.r, &tmp.c, &tmp.val);
    tmp.r -= 1;
    tmp.c -= 1;
    records[count] = tmp;
    ++count;
  }
  in.close();
}

void records_reorder_by_rows(const int &nz, record_t *records) {
  std::sort(records, records + nz, cmp_record_key_row);
}

void records_reorder_by_columns(const int &nz, record_t *records) {
  std::sort(records, records + nz, cmp_record_key_column);
}

void records_reorder_by_rows_columns(const int &nz, record_t *records) {
  std::sort(records, records + nz, cmp_record_key_column);
}

// source vector x random generation
void rand_gen(const int &len, double *x) {
  static thread_local std::mt19937 seed;
  std::uniform_real_distribution<> dis(-RAND_MAX, RAND_MAX);
  for (int i = 0; i < len; ++i) {
    x[i] = dis(seed);
    // fprintf(stdout, "%lf ", x[i]);
  }
}

// check the correctness of optimizer output and naive result
bool check(const size_t &len, double *output, double *result) {
  bool pass = true;
  for (std::size_t i = 0; i < len; ++i) {
    if ((output[i] - result[i]) > ESP) {
      pass = false;
      printf("fail at index %lu, out: %lf expected: %lf \n", i, output[i],
             result[i]);
      // break;
    }
  }
  return pass;
}

void cvt2bcsr(const int &Rows, const int &nz, const int &r, const int &c,
              const char* file_name, const int &block_nz,
              const std::vector<std::pair<int, int>> &block_idx,
              double *b_values, int *b_col_idx, int *b_row_start) {
  // 1. get b_row_start and b_col_idx
  b_row_start[0] = 0;
  b_col_idx[0] = block_idx[0].second;
  int before = block_idx[0].first;
  int current;
  int b_row_start_count = 0;
  for (int i = 1; i < block_idx.size(); ++i) {
    before = block_idx[i - 1].first;
    current = block_idx[i].first;
    b_col_idx[i] = block_idx[i].second;
    if (current != before) {
      b_row_start[++b_row_start_count] = i;
    }
  }
  b_row_start[Rows / r] = block_nz;
  // get index hash table
  using idx_hash_t = std::map<std::pair<int, int>, int>;
  idx_hash_t fuck_hash;
  for (int i = 0; i<block_idx.size(); ++i) {
    fuck_hash.insert(std::make_pair(block_idx[i], i));
  }
  // 2. get b_values
  std::ifstream in;
  in.open(file_name);
  if (!in.is_open()) {
    fprintf(stderr, "Fail to open file: %s\n", file_name);
  }
  int _M, _N, _nz;
  std::string line;
  while (std::getline(in, line) && line != "") {
    // skip the comments
    if (line[0] == '%')
      continue;
    sscanf(line.c_str(), "%d %d %d", &_M, &_N, &_nz);
    break;
  }
  int t_row, t_col;
  double t_val;
  for (int i = 0; i < nz; ++i) {
    std::getline(in, line);
    sscanf(line.c_str(), "%d %d %lf", &t_row, &t_col, &t_val);
    t_row -= 1;
    t_col -= 1;
    int _row = t_row / r;
    int _col = t_col / c;
    int off_row = t_row % r;
    int off_col = t_col % c;
    /*for (int j = 0; j < block_idx.size(); ++j) {
      if (block_idx[j] == std::make_pair(_row, _col)) {
        b_values[j * r * c + off_row * c + off_col] = t_val;
        break;
      }
    }*/
    int j = fuck_hash[std::make_pair(_row, _col)];
    b_values[j * r * c + off_row * c + off_col] = t_val;
  }
  in.close();
}

void cvt2csr(const int &rows, const int &nz, record_t *records, double *nz_vals,
             int *column_index, int *row_start) {
  // TODO: fill with a bit map
  int *row_tmp = (int *)malloc(rows * sizeof(int));
  if (!row_tmp) {
    fprintf(stdout, "%s:%d, Fail to malloc tmp!\n", __FILE__, __LINE__);
    return;
  }
  std::fill(row_tmp, row_tmp + rows, 0);
  int i;
  for (i = 0; i < nz; ++i) {
    ++row_tmp[records[i].r];
  }
  row_start[0] = 0;
  for (i = 1; i < rows + 1; ++i) {
    row_start[i] = row_start[i - 1] + row_tmp[i - 1];
  }
  free(row_tmp);
  row_start[rows] = nz;
  int j, one_row_count;
  // sort and get column index and non-zero values
  for (i = 0; i < rows; ++i) {
    one_row_count = row_start[i + 1] - row_start[i];
    record_t *tmp = (record_t *)malloc(one_row_count * sizeof(record_t));
    for (j = 0; j < one_row_count; ++j) {
      tmp[j] = records[row_start[i] + j];
    }
    // sort by column
    std::sort(tmp, tmp + one_row_count, cmp_record_key_column);
    for (j = 0; j < one_row_count; ++j) {
      column_index[row_start[i] + j] = tmp[j].c;
      nz_vals[row_start[i] + j] = tmp[j].val;
    }
    free(tmp);
  }
}

void cvt2csc(const int &columns, const int &nz, record_t *records,
             double *nz_vals, int *row_index, int *column_start) {
  // TODO: fill with a bit map
  int *col_tmp = (int *)malloc(columns * sizeof(int));
  if (!col_tmp) {
    fprintf(stdout, "%s:%d, Fail to malloc tmp!\n", __FILE__, __LINE__);
    return;
  }
  std::fill(col_tmp, col_tmp + columns, 0);
  int i;
  for (i = 0; i < nz; ++i) {
    ++col_tmp[records[i].c];
  }
  column_start[0] = 0;
  column_start[1] = col_tmp[0] + column_start[0];
  for (i = 1; i < columns + 1; ++i) {
    column_start[i] = column_start[i - 1] + col_tmp[i - 1];
  }
  free(col_tmp);
  column_start[columns] = nz;

  int j;
  int one_column_count;
  // sort and get column index and non-zero values
  for (i = 0; i < columns; ++i) {
    one_column_count = column_start[i + 1] - column_start[i];
    record_t *tmp = (record_t *)malloc(one_column_count * sizeof(record_t));
    for (j = 0; j < one_column_count; ++j) {
      tmp[j] = records[column_start[i] + j];
    }
    // sort by row
    std::sort(tmp, tmp + one_column_count, cmp_record_key_row);
    for (j = 0; j < one_column_count; ++j) {
      row_index[column_start[i] + j] = tmp[j].r;
      nz_vals[column_start[i] + j] = tmp[j].val;
    }
    free(tmp);
  }
}

#if 0
// convert the input matrix into csr format
void cvt2csr(const int &rows, const int &columns, const int &nz, double *A,
             double *nz_vals, int *column_index, int *row_start) {
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
  row_start[rows] = nz;
}

// convert the input matrix into csc format
void cvt2csc(const int &rows, const int &columns, const int &nz, double *A,
             double *nz_vals, int *row_index, int *column_start) {
  int count = 0;
  double element;
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
#endif // 0

#endif // _UTILS_H_
