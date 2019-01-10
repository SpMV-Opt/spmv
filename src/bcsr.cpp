/*
 *   Register Blocking CSR implement.
 *
 */
#include "spmv_opt.h"

BCSR_FP bcsr_routines[12][12][1] = {{{bcsr_1x1},
                                     {bcsr_1x2},
                                     {bcsr_1x3},
                                     {bcsr_1x4},
                                     {bcsr_1x5},
                                     {bcsr_1x6},
                                     {bcsr_1x7},
                                     {bcsr_1x8},
                                     {bcsr_1x9},
                                     {bcsr_1x10},
                                     {bcsr_1x11},
                                     {bcsr_1x12}},
                                    {{bcsr_2x1},
                                     {bcsr_2x2},
                                     {bcsr_2x3},
                                     {bcsr_2x4},
                                     {bcsr_2x5},
                                     {bcsr_2x6},
                                     {bcsr_2x7},
                                     {bcsr_2x8},
                                     {bcsr_2x9},
                                     {bcsr_2x10},
                                     {bcsr_2x11},
                                     {bcsr_2x12}},
                                    {{bcsr_3x1},
                                     {bcsr_3x2},
                                     {bcsr_3x3},
                                     {bcsr_3x4},
                                     {bcsr_3x5},
                                     {bcsr_3x6},
                                     {bcsr_3x7},
                                     {bcsr_3x8},
                                     {bcsr_3x9},
                                     {bcsr_3x10},
                                     {bcsr_3x11},
                                     {bcsr_3x12}},
                                    {{bcsr_4x1},
                                     {bcsr_4x2},
                                     {bcsr_4x3},
                                     {bcsr_4x4},
                                     {bcsr_4x5},
                                     {bcsr_4x6},
                                     {bcsr_4x7},
                                     {bcsr_4x8},
                                     {bcsr_4x9},
                                     {bcsr_4x10},
                                     {bcsr_4x11},
                                     {bcsr_4x12}},
                                    {{bcsr_5x1},
                                     {bcsr_5x2},
                                     {bcsr_5x3},
                                     {bcsr_5x4},
                                     {bcsr_5x5},
                                     {bcsr_5x6},
                                     {bcsr_5x7},
                                     {bcsr_5x8},
                                     {bcsr_5x9},
                                     {bcsr_5x10},
                                     {bcsr_5x11},
                                     {bcsr_5x12}},
                                    {{bcsr_6x1},
                                     {bcsr_6x2},
                                     {bcsr_6x3},
                                     {bcsr_6x4},
                                     {bcsr_6x5},
                                     {bcsr_6x6},
                                     {bcsr_6x7},
                                     {bcsr_6x8},
                                     {bcsr_6x9},
                                     {bcsr_6x10},
                                     {bcsr_6x11},
                                     {bcsr_6x12}},
                                    {{bcsr_7x1},
                                     {bcsr_7x2},
                                     {bcsr_7x3},
                                     {bcsr_7x4},
                                     {bcsr_7x5},
                                     {bcsr_7x6},
                                     {bcsr_7x7},
                                     {bcsr_7x8},
                                     {bcsr_7x9},
                                     {bcsr_7x10},
                                     {bcsr_7x11},
                                     {bcsr_7x12}},
                                    {{bcsr_8x1},
                                     {bcsr_8x2},
                                     {bcsr_8x3},
                                     {bcsr_8x4},
                                     {bcsr_8x5},
                                     {bcsr_8x6},
                                     {bcsr_8x7},
                                     {bcsr_8x8},
                                     {bcsr_8x9},
                                     {bcsr_8x10},
                                     {bcsr_8x11},
                                     {bcsr_8x12}},
                                    {{bcsr_9x1},
                                     {bcsr_9x2},
                                     {bcsr_9x3},
                                     {bcsr_9x4},
                                     {bcsr_9x5},
                                     {bcsr_9x6},
                                     {bcsr_9x7},
                                     {bcsr_9x8},
                                     {bcsr_9x9},
                                     {bcsr_9x10},
                                     {bcsr_9x11},
                                     {bcsr_9x12}},
                                    {{bcsr_10x1},
                                     {bcsr_10x2},
                                     {bcsr_10x3},
                                     {bcsr_10x4},
                                     {bcsr_10x5},
                                     {bcsr_10x6},
                                     {bcsr_10x7},
                                     {bcsr_10x8},
                                     {bcsr_10x9},
                                     {bcsr_10x10},
                                     {bcsr_10x11},
                                     {bcsr_10x12}},
                                    {{bcsr_11x1},
                                     {bcsr_11x2},
                                     {bcsr_11x3},
                                     {bcsr_11x4},
                                     {bcsr_11x5},
                                     {bcsr_11x6},
                                     {bcsr_11x7},
                                     {bcsr_11x8},
                                     {bcsr_11x9},
                                     {bcsr_11x10},
                                     {bcsr_11x11},
                                     {bcsr_11x12}},
                                    {{bcsr_12x1},
                                     {bcsr_12x2},
                                     {bcsr_12x3},
                                     {bcsr_12x4},
                                     {bcsr_12x5},
                                     {bcsr_12x6},
                                     {bcsr_12x7},
                                     {bcsr_12x8},
                                     {bcsr_12x9},
                                     {bcsr_12x10},
                                     {bcsr_12x11},
                                     {bcsr_12x12}}};

void bcsr_1x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      d0 += b_values[7] * x7;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      d0 += b_values[7] * x7;
      d0 += b_values[8] * x8;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      d0 += b_values[7] * x7;
      d0 += b_values[8] * x8;
      d0 += b_values[9] * x9;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      d0 += b_values[7] * x7;
      d0 += b_values[8] * x8;
      d0 += b_values[9] * x9;
      d0 += b_values[10] * x10;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_1x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[1 * i + 0];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 1 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d0 += b_values[1] * x1;
      d0 += b_values[2] * x2;
      d0 += b_values[3] * x3;
      d0 += b_values[4] * x4;
      d0 += b_values[5] * x5;
      d0 += b_values[6] * x6;
      d0 += b_values[7] * x7;
      d0 += b_values[8] * x8;
      d0 += b_values[9] * x9;
      d0 += b_values[10] * x10;
      d0 += b_values[11] * x11;
      y[1 * i + 0] = d0;
    }
  }
}

void bcsr_2x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_2x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[2 * i + 0];
    d1 = y[2 * i + 1];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 2 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      y[2 * i + 0] = d0;
      y[2 * i + 1] = d1;
    }
  }
}

void bcsr_3x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_3x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[3 * i + 0];
    d1 = y[3 * i + 1];
    d2 = y[3 * i + 2];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 3 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      y[3 * i + 0] = d0;
      y[3 * i + 1] = d1;
      y[3 * i + 2] = d2;
    }
  }
}

void bcsr_4x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_4x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[4 * i + 0];
    d1 = y[4 * i + 1];
    d2 = y[4 * i + 2];
    d3 = y[4 * i + 3];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 4 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      y[4 * i + 0] = d0;
      y[4 * i + 1] = d1;
      y[4 * i + 2] = d2;
      y[4 * i + 3] = d3;
    }
  }
}

void bcsr_5x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_5x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[5 * i + 0];
    d1 = y[5 * i + 1];
    d2 = y[5 * i + 2];
    d3 = y[5 * i + 3];
    d4 = y[5 * i + 4];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 5 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      y[5 * i + 0] = d0;
      y[5 * i + 1] = d1;
      y[5 * i + 2] = d2;
      y[5 * i + 3] = d3;
      y[5 * i + 4] = d4;
    }
  }
}

void bcsr_6x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_6x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
      x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[6 * i + 0];
    d1 = y[6 * i + 1];
    d2 = y[6 * i + 2];
    d3 = y[6 * i + 3];
    d4 = y[6 * i + 4];
    d5 = y[6 * i + 5];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 6 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      y[6 * i + 0] = d0;
      y[6 * i + 1] = d1;
      y[6 * i + 2] = d2;
      y[6 * i + 3] = d3;
      y[6 * i + 4] = d4;
      y[6 * i + 5] = d5;
    }
  }
}

void bcsr_7x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9,
      x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_7x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9,
      x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[7 * i + 0];
    d1 = y[7 * i + 1];
    d2 = y[7 * i + 2];
    d3 = y[7 * i + 3];
    d4 = y[7 * i + 4];
    d5 = y[7 * i + 5];
    d6 = y[7 * i + 6];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 7 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      y[7 * i + 0] = d0;
      y[7 * i + 1] = d1;
      y[7 * i + 2] = d2;
      y[7 * i + 3] = d3;
      y[7 * i + 4] = d4;
      y[7 * i + 5] = d5;
      y[7 * i + 6] = d6;
    }
  }
}

void bcsr_8x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      d7 += b_values[7] * x0;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d7 += b_values[14] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      d7 += b_values[15] * x1;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d7 += b_values[21] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d7 += b_values[22] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      d7 += b_values[23] * x2;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d7 += b_values[28] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d7 += b_values[29] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d7 += b_values[30] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      d7 += b_values[31] * x3;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d7 += b_values[35] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d7 += b_values[36] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d7 += b_values[37] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d7 += b_values[38] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      d7 += b_values[39] * x4;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d7 += b_values[42] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d7 += b_values[43] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d7 += b_values[44] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d7 += b_values[45] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d7 += b_values[46] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      d7 += b_values[47] * x5;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d7 += b_values[49] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d7 += b_values[50] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d7 += b_values[51] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d7 += b_values[52] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d7 += b_values[53] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d7 += b_values[54] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      d7 += b_values[55] * x6;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d7 += b_values[56] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d7 += b_values[57] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d7 += b_values[58] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d7 += b_values[59] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d7 += b_values[60] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d7 += b_values[61] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d7 += b_values[62] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      d7 += b_values[63] * x7;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d7 += b_values[63] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d7 += b_values[64] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d7 += b_values[65] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d7 += b_values[66] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d7 += b_values[67] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d7 += b_values[68] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d7 += b_values[69] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d7 += b_values[70] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      d7 += b_values[71] * x8;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d7 += b_values[70] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d7 += b_values[71] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d7 += b_values[72] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d7 += b_values[73] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d7 += b_values[74] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d7 += b_values[75] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d7 += b_values[76] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d7 += b_values[77] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d7 += b_values[78] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      d7 += b_values[79] * x9;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9,
      x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d7 += b_values[77] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d7 += b_values[78] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d7 += b_values[79] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d7 += b_values[80] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d7 += b_values[81] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d7 += b_values[82] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d7 += b_values[83] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d7 += b_values[84] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d7 += b_values[85] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d7 += b_values[86] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      d7 += b_values[87] * x10;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_8x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9,
      x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[8 * i + 0];
    d1 = y[8 * i + 1];
    d2 = y[8 * i + 2];
    d3 = y[8 * i + 3];
    d4 = y[8 * i + 4];
    d5 = y[8 * i + 5];
    d6 = y[8 * i + 6];
    d7 = y[8 * i + 7];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 8 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d7 += b_values[84] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d7 += b_values[85] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d7 += b_values[86] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d7 += b_values[87] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d7 += b_values[88] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d7 += b_values[89] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d7 += b_values[90] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d7 += b_values[91] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d7 += b_values[92] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d7 += b_values[93] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d7 += b_values[94] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      d7 += b_values[95] * x11;
      y[8 * i + 0] = d0;
      y[8 * i + 1] = d1;
      y[8 * i + 2] = d2;
      y[8 * i + 3] = d3;
      y[8 * i + 4] = d4;
      y[8 * i + 5] = d5;
      y[8 * i + 6] = d6;
      y[8 * i + 7] = d7;
    }
  }
}

void bcsr_9x1(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      d7 += b_values[7] * x0;
      d8 += b_values[8] * x0;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x2(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d7 += b_values[14] * x0;
      d8 += b_values[16] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      d7 += b_values[15] * x1;
      d8 += b_values[17] * x1;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x3(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d7 += b_values[21] * x0;
      d8 += b_values[24] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d7 += b_values[22] * x1;
      d8 += b_values[25] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      d7 += b_values[23] * x2;
      d8 += b_values[26] * x2;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x4(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d7 += b_values[28] * x0;
      d8 += b_values[32] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d7 += b_values[29] * x1;
      d8 += b_values[33] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d7 += b_values[30] * x2;
      d8 += b_values[34] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      d7 += b_values[31] * x3;
      d8 += b_values[35] * x3;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x5(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d7 += b_values[35] * x0;
      d8 += b_values[40] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d7 += b_values[36] * x1;
      d8 += b_values[41] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d7 += b_values[37] * x2;
      d8 += b_values[42] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d7 += b_values[38] * x3;
      d8 += b_values[43] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      d7 += b_values[39] * x4;
      d8 += b_values[44] * x4;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x6(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d7 += b_values[42] * x0;
      d8 += b_values[48] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d7 += b_values[43] * x1;
      d8 += b_values[49] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d7 += b_values[44] * x2;
      d8 += b_values[50] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d7 += b_values[45] * x3;
      d8 += b_values[51] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d7 += b_values[46] * x4;
      d8 += b_values[52] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      d7 += b_values[47] * x5;
      d8 += b_values[53] * x5;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x7(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d7 += b_values[49] * x0;
      d8 += b_values[56] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d7 += b_values[50] * x1;
      d8 += b_values[57] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d7 += b_values[51] * x2;
      d8 += b_values[58] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d7 += b_values[52] * x3;
      d8 += b_values[59] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d7 += b_values[53] * x4;
      d8 += b_values[60] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d7 += b_values[54] * x5;
      d8 += b_values[61] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      d7 += b_values[55] * x6;
      d8 += b_values[62] * x6;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x8(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d7 += b_values[56] * x0;
      d8 += b_values[64] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d7 += b_values[57] * x1;
      d8 += b_values[65] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d7 += b_values[58] * x2;
      d8 += b_values[66] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d7 += b_values[59] * x3;
      d8 += b_values[67] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d7 += b_values[60] * x4;
      d8 += b_values[68] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d7 += b_values[61] * x5;
      d8 += b_values[69] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d7 += b_values[62] * x6;
      d8 += b_values[70] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      d7 += b_values[63] * x7;
      d8 += b_values[71] * x7;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x9(const int &bm, const int *b_row_start, const int *b_col_idx,
              const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d7 += b_values[63] * x0;
      d8 += b_values[72] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d7 += b_values[64] * x1;
      d8 += b_values[73] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d7 += b_values[65] * x2;
      d8 += b_values[74] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d7 += b_values[66] * x3;
      d8 += b_values[75] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d7 += b_values[67] * x4;
      d8 += b_values[76] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d7 += b_values[68] * x5;
      d8 += b_values[77] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d7 += b_values[69] * x6;
      d8 += b_values[78] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d7 += b_values[70] * x7;
      d8 += b_values[79] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      d7 += b_values[71] * x8;
      d8 += b_values[80] * x8;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x10(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6, x7, x8,
      x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d7 += b_values[70] * x0;
      d8 += b_values[80] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d7 += b_values[71] * x1;
      d8 += b_values[81] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d7 += b_values[72] * x2;
      d8 += b_values[82] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d7 += b_values[73] * x3;
      d8 += b_values[83] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d7 += b_values[74] * x4;
      d8 += b_values[84] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d7 += b_values[75] * x5;
      d8 += b_values[85] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d7 += b_values[76] * x6;
      d8 += b_values[86] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d7 += b_values[77] * x7;
      d8 += b_values[87] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d7 += b_values[78] * x8;
      d8 += b_values[88] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      d7 += b_values[79] * x9;
      d8 += b_values[89] * x9;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x11(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6, x7, x8,
      x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d7 += b_values[77] * x0;
      d8 += b_values[88] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d7 += b_values[78] * x1;
      d8 += b_values[89] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d7 += b_values[79] * x2;
      d8 += b_values[90] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d7 += b_values[80] * x3;
      d8 += b_values[91] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d7 += b_values[81] * x4;
      d8 += b_values[92] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d7 += b_values[82] * x5;
      d8 += b_values[93] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d7 += b_values[83] * x6;
      d8 += b_values[94] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d7 += b_values[84] * x7;
      d8 += b_values[95] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d7 += b_values[85] * x8;
      d8 += b_values[96] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d7 += b_values[86] * x9;
      d8 += b_values[97] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      d7 += b_values[87] * x10;
      d8 += b_values[98] * x10;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_9x12(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, x0, x1, x2, x3, x4, x5, x6, x7, x8,
      x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[9 * i + 0];
    d1 = y[9 * i + 1];
    d2 = y[9 * i + 2];
    d3 = y[9 * i + 3];
    d4 = y[9 * i + 4];
    d5 = y[9 * i + 5];
    d6 = y[9 * i + 6];
    d7 = y[9 * i + 7];
    d8 = y[9 * i + 8];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 9 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d7 += b_values[84] * x0;
      d8 += b_values[96] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d7 += b_values[85] * x1;
      d8 += b_values[97] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d7 += b_values[86] * x2;
      d8 += b_values[98] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d7 += b_values[87] * x3;
      d8 += b_values[99] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d7 += b_values[88] * x4;
      d8 += b_values[100] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d7 += b_values[89] * x5;
      d8 += b_values[101] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d7 += b_values[90] * x6;
      d8 += b_values[102] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d7 += b_values[91] * x7;
      d8 += b_values[103] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d7 += b_values[92] * x8;
      d8 += b_values[104] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d7 += b_values[93] * x9;
      d8 += b_values[105] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d7 += b_values[94] * x10;
      d8 += b_values[106] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      d7 += b_values[95] * x11;
      d8 += b_values[107] * x11;
      y[9 * i + 0] = d0;
      y[9 * i + 1] = d1;
      y[9 * i + 2] = d2;
      y[9 * i + 3] = d3;
      y[9 * i + 4] = d4;
      y[9 * i + 5] = d5;
      y[9 * i + 6] = d6;
      y[9 * i + 7] = d7;
      y[9 * i + 8] = d8;
    }
  }
}

void bcsr_10x1(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      d7 += b_values[7] * x0;
      d8 += b_values[8] * x0;
      d9 += b_values[9] * x0;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x2(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d7 += b_values[14] * x0;
      d8 += b_values[16] * x0;
      d9 += b_values[18] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      d7 += b_values[15] * x1;
      d8 += b_values[17] * x1;
      d9 += b_values[19] * x1;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x3(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d7 += b_values[21] * x0;
      d8 += b_values[24] * x0;
      d9 += b_values[27] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d7 += b_values[22] * x1;
      d8 += b_values[25] * x1;
      d9 += b_values[28] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      d7 += b_values[23] * x2;
      d8 += b_values[26] * x2;
      d9 += b_values[29] * x2;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x4(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d7 += b_values[28] * x0;
      d8 += b_values[32] * x0;
      d9 += b_values[36] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d7 += b_values[29] * x1;
      d8 += b_values[33] * x1;
      d9 += b_values[37] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d7 += b_values[30] * x2;
      d8 += b_values[34] * x2;
      d9 += b_values[38] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      d7 += b_values[31] * x3;
      d8 += b_values[35] * x3;
      d9 += b_values[39] * x3;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x5(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d7 += b_values[35] * x0;
      d8 += b_values[40] * x0;
      d9 += b_values[45] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d7 += b_values[36] * x1;
      d8 += b_values[41] * x1;
      d9 += b_values[46] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d7 += b_values[37] * x2;
      d8 += b_values[42] * x2;
      d9 += b_values[47] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d7 += b_values[38] * x3;
      d8 += b_values[43] * x3;
      d9 += b_values[48] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      d7 += b_values[39] * x4;
      d8 += b_values[44] * x4;
      d9 += b_values[49] * x4;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x6(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d7 += b_values[42] * x0;
      d8 += b_values[48] * x0;
      d9 += b_values[54] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d7 += b_values[43] * x1;
      d8 += b_values[49] * x1;
      d9 += b_values[55] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d7 += b_values[44] * x2;
      d8 += b_values[50] * x2;
      d9 += b_values[56] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d7 += b_values[45] * x3;
      d8 += b_values[51] * x3;
      d9 += b_values[57] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d7 += b_values[46] * x4;
      d8 += b_values[52] * x4;
      d9 += b_values[58] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      d7 += b_values[47] * x5;
      d8 += b_values[53] * x5;
      d9 += b_values[59] * x5;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x7(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d7 += b_values[49] * x0;
      d8 += b_values[56] * x0;
      d9 += b_values[63] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d7 += b_values[50] * x1;
      d8 += b_values[57] * x1;
      d9 += b_values[64] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d7 += b_values[51] * x2;
      d8 += b_values[58] * x2;
      d9 += b_values[65] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d7 += b_values[52] * x3;
      d8 += b_values[59] * x3;
      d9 += b_values[66] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d7 += b_values[53] * x4;
      d8 += b_values[60] * x4;
      d9 += b_values[67] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d7 += b_values[54] * x5;
      d8 += b_values[61] * x5;
      d9 += b_values[68] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      d7 += b_values[55] * x6;
      d8 += b_values[62] * x6;
      d9 += b_values[69] * x6;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x8(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d7 += b_values[56] * x0;
      d8 += b_values[64] * x0;
      d9 += b_values[72] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d7 += b_values[57] * x1;
      d8 += b_values[65] * x1;
      d9 += b_values[73] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d7 += b_values[58] * x2;
      d8 += b_values[66] * x2;
      d9 += b_values[74] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d7 += b_values[59] * x3;
      d8 += b_values[67] * x3;
      d9 += b_values[75] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d7 += b_values[60] * x4;
      d8 += b_values[68] * x4;
      d9 += b_values[76] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d7 += b_values[61] * x5;
      d8 += b_values[69] * x5;
      d9 += b_values[77] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d7 += b_values[62] * x6;
      d8 += b_values[70] * x6;
      d9 += b_values[78] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      d7 += b_values[63] * x7;
      d8 += b_values[71] * x7;
      d9 += b_values[79] * x7;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x9(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6, x7,
      x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d7 += b_values[63] * x0;
      d8 += b_values[72] * x0;
      d9 += b_values[81] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d7 += b_values[64] * x1;
      d8 += b_values[73] * x1;
      d9 += b_values[82] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d7 += b_values[65] * x2;
      d8 += b_values[74] * x2;
      d9 += b_values[83] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d7 += b_values[66] * x3;
      d8 += b_values[75] * x3;
      d9 += b_values[84] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d7 += b_values[67] * x4;
      d8 += b_values[76] * x4;
      d9 += b_values[85] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d7 += b_values[68] * x5;
      d8 += b_values[77] * x5;
      d9 += b_values[86] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d7 += b_values[69] * x6;
      d8 += b_values[78] * x6;
      d9 += b_values[87] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d7 += b_values[70] * x7;
      d8 += b_values[79] * x7;
      d9 += b_values[88] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      d7 += b_values[71] * x8;
      d8 += b_values[80] * x8;
      d9 += b_values[89] * x8;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x10(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6, x7,
      x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d7 += b_values[70] * x0;
      d8 += b_values[80] * x0;
      d9 += b_values[90] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d7 += b_values[71] * x1;
      d8 += b_values[81] * x1;
      d9 += b_values[91] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d7 += b_values[72] * x2;
      d8 += b_values[82] * x2;
      d9 += b_values[92] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d7 += b_values[73] * x3;
      d8 += b_values[83] * x3;
      d9 += b_values[93] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d7 += b_values[74] * x4;
      d8 += b_values[84] * x4;
      d9 += b_values[94] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d7 += b_values[75] * x5;
      d8 += b_values[85] * x5;
      d9 += b_values[95] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d7 += b_values[76] * x6;
      d8 += b_values[86] * x6;
      d9 += b_values[96] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d7 += b_values[77] * x7;
      d8 += b_values[87] * x7;
      d9 += b_values[97] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d7 += b_values[78] * x8;
      d8 += b_values[88] * x8;
      d9 += b_values[98] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      d7 += b_values[79] * x9;
      d8 += b_values[89] * x9;
      d9 += b_values[99] * x9;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x11(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6, x7,
      x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d7 += b_values[77] * x0;
      d8 += b_values[88] * x0;
      d9 += b_values[99] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d7 += b_values[78] * x1;
      d8 += b_values[89] * x1;
      d9 += b_values[100] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d7 += b_values[79] * x2;
      d8 += b_values[90] * x2;
      d9 += b_values[101] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d7 += b_values[80] * x3;
      d8 += b_values[91] * x3;
      d9 += b_values[102] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d7 += b_values[81] * x4;
      d8 += b_values[92] * x4;
      d9 += b_values[103] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d7 += b_values[82] * x5;
      d8 += b_values[93] * x5;
      d9 += b_values[104] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d7 += b_values[83] * x6;
      d8 += b_values[94] * x6;
      d9 += b_values[105] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d7 += b_values[84] * x7;
      d8 += b_values[95] * x7;
      d9 += b_values[106] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d7 += b_values[85] * x8;
      d8 += b_values[96] * x8;
      d9 += b_values[107] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d7 += b_values[86] * x9;
      d8 += b_values[97] * x9;
      d9 += b_values[108] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      d7 += b_values[87] * x10;
      d8 += b_values[98] * x10;
      d9 += b_values[109] * x10;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_10x12(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, x0, x1, x2, x3, x4, x5, x6, x7,
      x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[10 * i + 0];
    d1 = y[10 * i + 1];
    d2 = y[10 * i + 2];
    d3 = y[10 * i + 3];
    d4 = y[10 * i + 4];
    d5 = y[10 * i + 5];
    d6 = y[10 * i + 6];
    d7 = y[10 * i + 7];
    d8 = y[10 * i + 8];
    d9 = y[10 * i + 9];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 10 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d7 += b_values[84] * x0;
      d8 += b_values[96] * x0;
      d9 += b_values[108] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d7 += b_values[85] * x1;
      d8 += b_values[97] * x1;
      d9 += b_values[109] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d7 += b_values[86] * x2;
      d8 += b_values[98] * x2;
      d9 += b_values[110] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d7 += b_values[87] * x3;
      d8 += b_values[99] * x3;
      d9 += b_values[111] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d7 += b_values[88] * x4;
      d8 += b_values[100] * x4;
      d9 += b_values[112] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d7 += b_values[89] * x5;
      d8 += b_values[101] * x5;
      d9 += b_values[113] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d7 += b_values[90] * x6;
      d8 += b_values[102] * x6;
      d9 += b_values[114] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d7 += b_values[91] * x7;
      d8 += b_values[103] * x7;
      d9 += b_values[115] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d7 += b_values[92] * x8;
      d8 += b_values[104] * x8;
      d9 += b_values[116] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d7 += b_values[93] * x9;
      d8 += b_values[105] * x9;
      d9 += b_values[117] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d7 += b_values[94] * x10;
      d8 += b_values[106] * x10;
      d9 += b_values[118] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      d7 += b_values[95] * x11;
      d8 += b_values[107] * x11;
      d9 += b_values[119] * x11;
      y[10 * i + 0] = d0;
      y[10 * i + 1] = d1;
      y[10 * i + 2] = d2;
      y[10 * i + 3] = d3;
      y[10 * i + 4] = d4;
      y[10 * i + 5] = d5;
      y[10 * i + 6] = d6;
      y[10 * i + 7] = d7;
      y[10 * i + 8] = d8;
      y[10 * i + 9] = d9;
    }
  }
}

void bcsr_11x1(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      d7 += b_values[7] * x0;
      d8 += b_values[8] * x0;
      d9 += b_values[9] * x0;
      d10 += b_values[10] * x0;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x2(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d7 += b_values[14] * x0;
      d8 += b_values[16] * x0;
      d9 += b_values[18] * x0;
      d10 += b_values[20] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      d7 += b_values[15] * x1;
      d8 += b_values[17] * x1;
      d9 += b_values[19] * x1;
      d10 += b_values[21] * x1;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x3(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d7 += b_values[21] * x0;
      d8 += b_values[24] * x0;
      d9 += b_values[27] * x0;
      d10 += b_values[30] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d7 += b_values[22] * x1;
      d8 += b_values[25] * x1;
      d9 += b_values[28] * x1;
      d10 += b_values[31] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      d7 += b_values[23] * x2;
      d8 += b_values[26] * x2;
      d9 += b_values[29] * x2;
      d10 += b_values[32] * x2;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x4(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d7 += b_values[28] * x0;
      d8 += b_values[32] * x0;
      d9 += b_values[36] * x0;
      d10 += b_values[40] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d7 += b_values[29] * x1;
      d8 += b_values[33] * x1;
      d9 += b_values[37] * x1;
      d10 += b_values[41] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d7 += b_values[30] * x2;
      d8 += b_values[34] * x2;
      d9 += b_values[38] * x2;
      d10 += b_values[42] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      d7 += b_values[31] * x3;
      d8 += b_values[35] * x3;
      d9 += b_values[39] * x3;
      d10 += b_values[43] * x3;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x5(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d7 += b_values[35] * x0;
      d8 += b_values[40] * x0;
      d9 += b_values[45] * x0;
      d10 += b_values[50] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d7 += b_values[36] * x1;
      d8 += b_values[41] * x1;
      d9 += b_values[46] * x1;
      d10 += b_values[51] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d7 += b_values[37] * x2;
      d8 += b_values[42] * x2;
      d9 += b_values[47] * x2;
      d10 += b_values[52] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d7 += b_values[38] * x3;
      d8 += b_values[43] * x3;
      d9 += b_values[48] * x3;
      d10 += b_values[53] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      d7 += b_values[39] * x4;
      d8 += b_values[44] * x4;
      d9 += b_values[49] * x4;
      d10 += b_values[54] * x4;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x6(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d7 += b_values[42] * x0;
      d8 += b_values[48] * x0;
      d9 += b_values[54] * x0;
      d10 += b_values[60] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d7 += b_values[43] * x1;
      d8 += b_values[49] * x1;
      d9 += b_values[55] * x1;
      d10 += b_values[61] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d7 += b_values[44] * x2;
      d8 += b_values[50] * x2;
      d9 += b_values[56] * x2;
      d10 += b_values[62] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d7 += b_values[45] * x3;
      d8 += b_values[51] * x3;
      d9 += b_values[57] * x3;
      d10 += b_values[63] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d7 += b_values[46] * x4;
      d8 += b_values[52] * x4;
      d9 += b_values[58] * x4;
      d10 += b_values[64] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      d7 += b_values[47] * x5;
      d8 += b_values[53] * x5;
      d9 += b_values[59] * x5;
      d10 += b_values[65] * x5;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x7(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d7 += b_values[49] * x0;
      d8 += b_values[56] * x0;
      d9 += b_values[63] * x0;
      d10 += b_values[70] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d7 += b_values[50] * x1;
      d8 += b_values[57] * x1;
      d9 += b_values[64] * x1;
      d10 += b_values[71] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d7 += b_values[51] * x2;
      d8 += b_values[58] * x2;
      d9 += b_values[65] * x2;
      d10 += b_values[72] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d7 += b_values[52] * x3;
      d8 += b_values[59] * x3;
      d9 += b_values[66] * x3;
      d10 += b_values[73] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d7 += b_values[53] * x4;
      d8 += b_values[60] * x4;
      d9 += b_values[67] * x4;
      d10 += b_values[74] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d7 += b_values[54] * x5;
      d8 += b_values[61] * x5;
      d9 += b_values[68] * x5;
      d10 += b_values[75] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      d7 += b_values[55] * x6;
      d8 += b_values[62] * x6;
      d9 += b_values[69] * x6;
      d10 += b_values[76] * x6;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x8(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d7 += b_values[56] * x0;
      d8 += b_values[64] * x0;
      d9 += b_values[72] * x0;
      d10 += b_values[80] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d7 += b_values[57] * x1;
      d8 += b_values[65] * x1;
      d9 += b_values[73] * x1;
      d10 += b_values[81] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d7 += b_values[58] * x2;
      d8 += b_values[66] * x2;
      d9 += b_values[74] * x2;
      d10 += b_values[82] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d7 += b_values[59] * x3;
      d8 += b_values[67] * x3;
      d9 += b_values[75] * x3;
      d10 += b_values[83] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d7 += b_values[60] * x4;
      d8 += b_values[68] * x4;
      d9 += b_values[76] * x4;
      d10 += b_values[84] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d7 += b_values[61] * x5;
      d8 += b_values[69] * x5;
      d9 += b_values[77] * x5;
      d10 += b_values[85] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d7 += b_values[62] * x6;
      d8 += b_values[70] * x6;
      d9 += b_values[78] * x6;
      d10 += b_values[86] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      d7 += b_values[63] * x7;
      d8 += b_values[71] * x7;
      d9 += b_values[79] * x7;
      d10 += b_values[87] * x7;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x9(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d7 += b_values[63] * x0;
      d8 += b_values[72] * x0;
      d9 += b_values[81] * x0;
      d10 += b_values[90] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d7 += b_values[64] * x1;
      d8 += b_values[73] * x1;
      d9 += b_values[82] * x1;
      d10 += b_values[91] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d7 += b_values[65] * x2;
      d8 += b_values[74] * x2;
      d9 += b_values[83] * x2;
      d10 += b_values[92] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d7 += b_values[66] * x3;
      d8 += b_values[75] * x3;
      d9 += b_values[84] * x3;
      d10 += b_values[93] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d7 += b_values[67] * x4;
      d8 += b_values[76] * x4;
      d9 += b_values[85] * x4;
      d10 += b_values[94] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d7 += b_values[68] * x5;
      d8 += b_values[77] * x5;
      d9 += b_values[86] * x5;
      d10 += b_values[95] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d7 += b_values[69] * x6;
      d8 += b_values[78] * x6;
      d9 += b_values[87] * x6;
      d10 += b_values[96] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d7 += b_values[70] * x7;
      d8 += b_values[79] * x7;
      d9 += b_values[88] * x7;
      d10 += b_values[97] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      d7 += b_values[71] * x8;
      d8 += b_values[80] * x8;
      d9 += b_values[89] * x8;
      d10 += b_values[98] * x8;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x10(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d7 += b_values[70] * x0;
      d8 += b_values[80] * x0;
      d9 += b_values[90] * x0;
      d10 += b_values[100] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d7 += b_values[71] * x1;
      d8 += b_values[81] * x1;
      d9 += b_values[91] * x1;
      d10 += b_values[101] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d7 += b_values[72] * x2;
      d8 += b_values[82] * x2;
      d9 += b_values[92] * x2;
      d10 += b_values[102] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d7 += b_values[73] * x3;
      d8 += b_values[83] * x3;
      d9 += b_values[93] * x3;
      d10 += b_values[103] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d7 += b_values[74] * x4;
      d8 += b_values[84] * x4;
      d9 += b_values[94] * x4;
      d10 += b_values[104] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d7 += b_values[75] * x5;
      d8 += b_values[85] * x5;
      d9 += b_values[95] * x5;
      d10 += b_values[105] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d7 += b_values[76] * x6;
      d8 += b_values[86] * x6;
      d9 += b_values[96] * x6;
      d10 += b_values[106] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d7 += b_values[77] * x7;
      d8 += b_values[87] * x7;
      d9 += b_values[97] * x7;
      d10 += b_values[107] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d7 += b_values[78] * x8;
      d8 += b_values[88] * x8;
      d9 += b_values[98] * x8;
      d10 += b_values[108] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      d7 += b_values[79] * x9;
      d8 += b_values[89] * x9;
      d9 += b_values[99] * x9;
      d10 += b_values[109] * x9;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x11(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d7 += b_values[77] * x0;
      d8 += b_values[88] * x0;
      d9 += b_values[99] * x0;
      d10 += b_values[110] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d7 += b_values[78] * x1;
      d8 += b_values[89] * x1;
      d9 += b_values[100] * x1;
      d10 += b_values[111] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d7 += b_values[79] * x2;
      d8 += b_values[90] * x2;
      d9 += b_values[101] * x2;
      d10 += b_values[112] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d7 += b_values[80] * x3;
      d8 += b_values[91] * x3;
      d9 += b_values[102] * x3;
      d10 += b_values[113] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d7 += b_values[81] * x4;
      d8 += b_values[92] * x4;
      d9 += b_values[103] * x4;
      d10 += b_values[114] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d7 += b_values[82] * x5;
      d8 += b_values[93] * x5;
      d9 += b_values[104] * x5;
      d10 += b_values[115] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d7 += b_values[83] * x6;
      d8 += b_values[94] * x6;
      d9 += b_values[105] * x6;
      d10 += b_values[116] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d7 += b_values[84] * x7;
      d8 += b_values[95] * x7;
      d9 += b_values[106] * x7;
      d10 += b_values[117] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d7 += b_values[85] * x8;
      d8 += b_values[96] * x8;
      d9 += b_values[107] * x8;
      d10 += b_values[118] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d7 += b_values[86] * x9;
      d8 += b_values[97] * x9;
      d9 += b_values[108] * x9;
      d10 += b_values[119] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      d7 += b_values[87] * x10;
      d8 += b_values[98] * x10;
      d9 += b_values[109] * x10;
      d10 += b_values[120] * x10;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_11x12(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, x0, x1, x2, x3, x4, x5,
      x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[11 * i + 0];
    d1 = y[11 * i + 1];
    d2 = y[11 * i + 2];
    d3 = y[11 * i + 3];
    d4 = y[11 * i + 4];
    d5 = y[11 * i + 5];
    d6 = y[11 * i + 6];
    d7 = y[11 * i + 7];
    d8 = y[11 * i + 8];
    d9 = y[11 * i + 9];
    d10 = y[11 * i + 10];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 11 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d7 += b_values[84] * x0;
      d8 += b_values[96] * x0;
      d9 += b_values[108] * x0;
      d10 += b_values[120] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d7 += b_values[85] * x1;
      d8 += b_values[97] * x1;
      d9 += b_values[109] * x1;
      d10 += b_values[121] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d7 += b_values[86] * x2;
      d8 += b_values[98] * x2;
      d9 += b_values[110] * x2;
      d10 += b_values[122] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d7 += b_values[87] * x3;
      d8 += b_values[99] * x3;
      d9 += b_values[111] * x3;
      d10 += b_values[123] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d7 += b_values[88] * x4;
      d8 += b_values[100] * x4;
      d9 += b_values[112] * x4;
      d10 += b_values[124] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d7 += b_values[89] * x5;
      d8 += b_values[101] * x5;
      d9 += b_values[113] * x5;
      d10 += b_values[125] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d7 += b_values[90] * x6;
      d8 += b_values[102] * x6;
      d9 += b_values[114] * x6;
      d10 += b_values[126] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d7 += b_values[91] * x7;
      d8 += b_values[103] * x7;
      d9 += b_values[115] * x7;
      d10 += b_values[127] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d7 += b_values[92] * x8;
      d8 += b_values[104] * x8;
      d9 += b_values[116] * x8;
      d10 += b_values[128] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d7 += b_values[93] * x9;
      d8 += b_values[105] * x9;
      d9 += b_values[117] * x9;
      d10 += b_values[129] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d7 += b_values[94] * x10;
      d8 += b_values[106] * x10;
      d9 += b_values[118] * x10;
      d10 += b_values[130] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      d7 += b_values[95] * x11;
      d8 += b_values[107] * x11;
      d9 += b_values[119] * x11;
      d10 += b_values[131] * x11;
      y[11 * i + 0] = d0;
      y[11 * i + 1] = d1;
      y[11 * i + 2] = d2;
      y[11 * i + 3] = d3;
      y[11 * i + 4] = d4;
      y[11 * i + 5] = d5;
      y[11 * i + 6] = d6;
      y[11 * i + 7] = d7;
      y[11 * i + 8] = d8;
      y[11 * i + 9] = d9;
      y[11 * i + 10] = d10;
    }
  }
}

void bcsr_12x1(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 1) {
      x0 = x[1 * b_col_idx[j] + 0];
      d0 += b_values[0] * x0;
      d1 += b_values[1] * x0;
      d2 += b_values[2] * x0;
      d3 += b_values[3] * x0;
      d4 += b_values[4] * x0;
      d5 += b_values[5] * x0;
      d6 += b_values[6] * x0;
      d7 += b_values[7] * x0;
      d8 += b_values[8] * x0;
      d9 += b_values[9] * x0;
      d10 += b_values[10] * x0;
      d11 += b_values[11] * x0;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x2(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 2) {
      x0 = x[2 * b_col_idx[j] + 0];
      x1 = x[2 * b_col_idx[j] + 1];
      d0 += b_values[0] * x0;
      d1 += b_values[2] * x0;
      d2 += b_values[4] * x0;
      d3 += b_values[6] * x0;
      d4 += b_values[8] * x0;
      d5 += b_values[10] * x0;
      d6 += b_values[12] * x0;
      d7 += b_values[14] * x0;
      d8 += b_values[16] * x0;
      d9 += b_values[18] * x0;
      d10 += b_values[20] * x0;
      d11 += b_values[22] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[3] * x1;
      d2 += b_values[5] * x1;
      d3 += b_values[7] * x1;
      d4 += b_values[9] * x1;
      d5 += b_values[11] * x1;
      d6 += b_values[13] * x1;
      d7 += b_values[15] * x1;
      d8 += b_values[17] * x1;
      d9 += b_values[19] * x1;
      d10 += b_values[21] * x1;
      d11 += b_values[23] * x1;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x3(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 3) {
      x0 = x[3 * b_col_idx[j] + 0];
      x1 = x[3 * b_col_idx[j] + 1];
      x2 = x[3 * b_col_idx[j] + 2];
      d0 += b_values[0] * x0;
      d1 += b_values[3] * x0;
      d2 += b_values[6] * x0;
      d3 += b_values[9] * x0;
      d4 += b_values[12] * x0;
      d5 += b_values[15] * x0;
      d6 += b_values[18] * x0;
      d7 += b_values[21] * x0;
      d8 += b_values[24] * x0;
      d9 += b_values[27] * x0;
      d10 += b_values[30] * x0;
      d11 += b_values[33] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[4] * x1;
      d2 += b_values[7] * x1;
      d3 += b_values[10] * x1;
      d4 += b_values[13] * x1;
      d5 += b_values[16] * x1;
      d6 += b_values[19] * x1;
      d7 += b_values[22] * x1;
      d8 += b_values[25] * x1;
      d9 += b_values[28] * x1;
      d10 += b_values[31] * x1;
      d11 += b_values[34] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[5] * x2;
      d2 += b_values[8] * x2;
      d3 += b_values[11] * x2;
      d4 += b_values[14] * x2;
      d5 += b_values[17] * x2;
      d6 += b_values[20] * x2;
      d7 += b_values[23] * x2;
      d8 += b_values[26] * x2;
      d9 += b_values[29] * x2;
      d10 += b_values[32] * x2;
      d11 += b_values[35] * x2;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x4(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 4) {
      x0 = x[4 * b_col_idx[j] + 0];
      x1 = x[4 * b_col_idx[j] + 1];
      x2 = x[4 * b_col_idx[j] + 2];
      x3 = x[4 * b_col_idx[j] + 3];
      d0 += b_values[0] * x0;
      d1 += b_values[4] * x0;
      d2 += b_values[8] * x0;
      d3 += b_values[12] * x0;
      d4 += b_values[16] * x0;
      d5 += b_values[20] * x0;
      d6 += b_values[24] * x0;
      d7 += b_values[28] * x0;
      d8 += b_values[32] * x0;
      d9 += b_values[36] * x0;
      d10 += b_values[40] * x0;
      d11 += b_values[44] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[5] * x1;
      d2 += b_values[9] * x1;
      d3 += b_values[13] * x1;
      d4 += b_values[17] * x1;
      d5 += b_values[21] * x1;
      d6 += b_values[25] * x1;
      d7 += b_values[29] * x1;
      d8 += b_values[33] * x1;
      d9 += b_values[37] * x1;
      d10 += b_values[41] * x1;
      d11 += b_values[45] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[6] * x2;
      d2 += b_values[10] * x2;
      d3 += b_values[14] * x2;
      d4 += b_values[18] * x2;
      d5 += b_values[22] * x2;
      d6 += b_values[26] * x2;
      d7 += b_values[30] * x2;
      d8 += b_values[34] * x2;
      d9 += b_values[38] * x2;
      d10 += b_values[42] * x2;
      d11 += b_values[46] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[7] * x3;
      d2 += b_values[11] * x3;
      d3 += b_values[15] * x3;
      d4 += b_values[19] * x3;
      d5 += b_values[23] * x3;
      d6 += b_values[27] * x3;
      d7 += b_values[31] * x3;
      d8 += b_values[35] * x3;
      d9 += b_values[39] * x3;
      d10 += b_values[43] * x3;
      d11 += b_values[47] * x3;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x5(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 5) {
      x0 = x[5 * b_col_idx[j] + 0];
      x1 = x[5 * b_col_idx[j] + 1];
      x2 = x[5 * b_col_idx[j] + 2];
      x3 = x[5 * b_col_idx[j] + 3];
      x4 = x[5 * b_col_idx[j] + 4];
      d0 += b_values[0] * x0;
      d1 += b_values[5] * x0;
      d2 += b_values[10] * x0;
      d3 += b_values[15] * x0;
      d4 += b_values[20] * x0;
      d5 += b_values[25] * x0;
      d6 += b_values[30] * x0;
      d7 += b_values[35] * x0;
      d8 += b_values[40] * x0;
      d9 += b_values[45] * x0;
      d10 += b_values[50] * x0;
      d11 += b_values[55] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[6] * x1;
      d2 += b_values[11] * x1;
      d3 += b_values[16] * x1;
      d4 += b_values[21] * x1;
      d5 += b_values[26] * x1;
      d6 += b_values[31] * x1;
      d7 += b_values[36] * x1;
      d8 += b_values[41] * x1;
      d9 += b_values[46] * x1;
      d10 += b_values[51] * x1;
      d11 += b_values[56] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[7] * x2;
      d2 += b_values[12] * x2;
      d3 += b_values[17] * x2;
      d4 += b_values[22] * x2;
      d5 += b_values[27] * x2;
      d6 += b_values[32] * x2;
      d7 += b_values[37] * x2;
      d8 += b_values[42] * x2;
      d9 += b_values[47] * x2;
      d10 += b_values[52] * x2;
      d11 += b_values[57] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[8] * x3;
      d2 += b_values[13] * x3;
      d3 += b_values[18] * x3;
      d4 += b_values[23] * x3;
      d5 += b_values[28] * x3;
      d6 += b_values[33] * x3;
      d7 += b_values[38] * x3;
      d8 += b_values[43] * x3;
      d9 += b_values[48] * x3;
      d10 += b_values[53] * x3;
      d11 += b_values[58] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[9] * x4;
      d2 += b_values[14] * x4;
      d3 += b_values[19] * x4;
      d4 += b_values[24] * x4;
      d5 += b_values[29] * x4;
      d6 += b_values[34] * x4;
      d7 += b_values[39] * x4;
      d8 += b_values[44] * x4;
      d9 += b_values[49] * x4;
      d10 += b_values[54] * x4;
      d11 += b_values[59] * x4;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x6(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 6) {
      x0 = x[6 * b_col_idx[j] + 0];
      x1 = x[6 * b_col_idx[j] + 1];
      x2 = x[6 * b_col_idx[j] + 2];
      x3 = x[6 * b_col_idx[j] + 3];
      x4 = x[6 * b_col_idx[j] + 4];
      x5 = x[6 * b_col_idx[j] + 5];
      d0 += b_values[0] * x0;
      d1 += b_values[6] * x0;
      d2 += b_values[12] * x0;
      d3 += b_values[18] * x0;
      d4 += b_values[24] * x0;
      d5 += b_values[30] * x0;
      d6 += b_values[36] * x0;
      d7 += b_values[42] * x0;
      d8 += b_values[48] * x0;
      d9 += b_values[54] * x0;
      d10 += b_values[60] * x0;
      d11 += b_values[66] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[7] * x1;
      d2 += b_values[13] * x1;
      d3 += b_values[19] * x1;
      d4 += b_values[25] * x1;
      d5 += b_values[31] * x1;
      d6 += b_values[37] * x1;
      d7 += b_values[43] * x1;
      d8 += b_values[49] * x1;
      d9 += b_values[55] * x1;
      d10 += b_values[61] * x1;
      d11 += b_values[67] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[8] * x2;
      d2 += b_values[14] * x2;
      d3 += b_values[20] * x2;
      d4 += b_values[26] * x2;
      d5 += b_values[32] * x2;
      d6 += b_values[38] * x2;
      d7 += b_values[44] * x2;
      d8 += b_values[50] * x2;
      d9 += b_values[56] * x2;
      d10 += b_values[62] * x2;
      d11 += b_values[68] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[9] * x3;
      d2 += b_values[15] * x3;
      d3 += b_values[21] * x3;
      d4 += b_values[27] * x3;
      d5 += b_values[33] * x3;
      d6 += b_values[39] * x3;
      d7 += b_values[45] * x3;
      d8 += b_values[51] * x3;
      d9 += b_values[57] * x3;
      d10 += b_values[63] * x3;
      d11 += b_values[69] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[10] * x4;
      d2 += b_values[16] * x4;
      d3 += b_values[22] * x4;
      d4 += b_values[28] * x4;
      d5 += b_values[34] * x4;
      d6 += b_values[40] * x4;
      d7 += b_values[46] * x4;
      d8 += b_values[52] * x4;
      d9 += b_values[58] * x4;
      d10 += b_values[64] * x4;
      d11 += b_values[70] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[11] * x5;
      d2 += b_values[17] * x5;
      d3 += b_values[23] * x5;
      d4 += b_values[29] * x5;
      d5 += b_values[35] * x5;
      d6 += b_values[41] * x5;
      d7 += b_values[47] * x5;
      d8 += b_values[53] * x5;
      d9 += b_values[59] * x5;
      d10 += b_values[65] * x5;
      d11 += b_values[71] * x5;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x7(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 7) {
      x0 = x[7 * b_col_idx[j] + 0];
      x1 = x[7 * b_col_idx[j] + 1];
      x2 = x[7 * b_col_idx[j] + 2];
      x3 = x[7 * b_col_idx[j] + 3];
      x4 = x[7 * b_col_idx[j] + 4];
      x5 = x[7 * b_col_idx[j] + 5];
      x6 = x[7 * b_col_idx[j] + 6];
      d0 += b_values[0] * x0;
      d1 += b_values[7] * x0;
      d2 += b_values[14] * x0;
      d3 += b_values[21] * x0;
      d4 += b_values[28] * x0;
      d5 += b_values[35] * x0;
      d6 += b_values[42] * x0;
      d7 += b_values[49] * x0;
      d8 += b_values[56] * x0;
      d9 += b_values[63] * x0;
      d10 += b_values[70] * x0;
      d11 += b_values[77] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[8] * x1;
      d2 += b_values[15] * x1;
      d3 += b_values[22] * x1;
      d4 += b_values[29] * x1;
      d5 += b_values[36] * x1;
      d6 += b_values[43] * x1;
      d7 += b_values[50] * x1;
      d8 += b_values[57] * x1;
      d9 += b_values[64] * x1;
      d10 += b_values[71] * x1;
      d11 += b_values[78] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[9] * x2;
      d2 += b_values[16] * x2;
      d3 += b_values[23] * x2;
      d4 += b_values[30] * x2;
      d5 += b_values[37] * x2;
      d6 += b_values[44] * x2;
      d7 += b_values[51] * x2;
      d8 += b_values[58] * x2;
      d9 += b_values[65] * x2;
      d10 += b_values[72] * x2;
      d11 += b_values[79] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[10] * x3;
      d2 += b_values[17] * x3;
      d3 += b_values[24] * x3;
      d4 += b_values[31] * x3;
      d5 += b_values[38] * x3;
      d6 += b_values[45] * x3;
      d7 += b_values[52] * x3;
      d8 += b_values[59] * x3;
      d9 += b_values[66] * x3;
      d10 += b_values[73] * x3;
      d11 += b_values[80] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[11] * x4;
      d2 += b_values[18] * x4;
      d3 += b_values[25] * x4;
      d4 += b_values[32] * x4;
      d5 += b_values[39] * x4;
      d6 += b_values[46] * x4;
      d7 += b_values[53] * x4;
      d8 += b_values[60] * x4;
      d9 += b_values[67] * x4;
      d10 += b_values[74] * x4;
      d11 += b_values[81] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[12] * x5;
      d2 += b_values[19] * x5;
      d3 += b_values[26] * x5;
      d4 += b_values[33] * x5;
      d5 += b_values[40] * x5;
      d6 += b_values[47] * x5;
      d7 += b_values[54] * x5;
      d8 += b_values[61] * x5;
      d9 += b_values[68] * x5;
      d10 += b_values[75] * x5;
      d11 += b_values[82] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[13] * x6;
      d2 += b_values[20] * x6;
      d3 += b_values[27] * x6;
      d4 += b_values[34] * x6;
      d5 += b_values[41] * x6;
      d6 += b_values[48] * x6;
      d7 += b_values[55] * x6;
      d8 += b_values[62] * x6;
      d9 += b_values[69] * x6;
      d10 += b_values[76] * x6;
      d11 += b_values[83] * x6;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x8(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6, x7;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 8) {
      x0 = x[8 * b_col_idx[j] + 0];
      x1 = x[8 * b_col_idx[j] + 1];
      x2 = x[8 * b_col_idx[j] + 2];
      x3 = x[8 * b_col_idx[j] + 3];
      x4 = x[8 * b_col_idx[j] + 4];
      x5 = x[8 * b_col_idx[j] + 5];
      x6 = x[8 * b_col_idx[j] + 6];
      x7 = x[8 * b_col_idx[j] + 7];
      d0 += b_values[0] * x0;
      d1 += b_values[8] * x0;
      d2 += b_values[16] * x0;
      d3 += b_values[24] * x0;
      d4 += b_values[32] * x0;
      d5 += b_values[40] * x0;
      d6 += b_values[48] * x0;
      d7 += b_values[56] * x0;
      d8 += b_values[64] * x0;
      d9 += b_values[72] * x0;
      d10 += b_values[80] * x0;
      d11 += b_values[88] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[9] * x1;
      d2 += b_values[17] * x1;
      d3 += b_values[25] * x1;
      d4 += b_values[33] * x1;
      d5 += b_values[41] * x1;
      d6 += b_values[49] * x1;
      d7 += b_values[57] * x1;
      d8 += b_values[65] * x1;
      d9 += b_values[73] * x1;
      d10 += b_values[81] * x1;
      d11 += b_values[89] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[10] * x2;
      d2 += b_values[18] * x2;
      d3 += b_values[26] * x2;
      d4 += b_values[34] * x2;
      d5 += b_values[42] * x2;
      d6 += b_values[50] * x2;
      d7 += b_values[58] * x2;
      d8 += b_values[66] * x2;
      d9 += b_values[74] * x2;
      d10 += b_values[82] * x2;
      d11 += b_values[90] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[11] * x3;
      d2 += b_values[19] * x3;
      d3 += b_values[27] * x3;
      d4 += b_values[35] * x3;
      d5 += b_values[43] * x3;
      d6 += b_values[51] * x3;
      d7 += b_values[59] * x3;
      d8 += b_values[67] * x3;
      d9 += b_values[75] * x3;
      d10 += b_values[83] * x3;
      d11 += b_values[91] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[12] * x4;
      d2 += b_values[20] * x4;
      d3 += b_values[28] * x4;
      d4 += b_values[36] * x4;
      d5 += b_values[44] * x4;
      d6 += b_values[52] * x4;
      d7 += b_values[60] * x4;
      d8 += b_values[68] * x4;
      d9 += b_values[76] * x4;
      d10 += b_values[84] * x4;
      d11 += b_values[92] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[13] * x5;
      d2 += b_values[21] * x5;
      d3 += b_values[29] * x5;
      d4 += b_values[37] * x5;
      d5 += b_values[45] * x5;
      d6 += b_values[53] * x5;
      d7 += b_values[61] * x5;
      d8 += b_values[69] * x5;
      d9 += b_values[77] * x5;
      d10 += b_values[85] * x5;
      d11 += b_values[93] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[14] * x6;
      d2 += b_values[22] * x6;
      d3 += b_values[30] * x6;
      d4 += b_values[38] * x6;
      d5 += b_values[46] * x6;
      d6 += b_values[54] * x6;
      d7 += b_values[62] * x6;
      d8 += b_values[70] * x6;
      d9 += b_values[78] * x6;
      d10 += b_values[86] * x6;
      d11 += b_values[94] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[15] * x7;
      d2 += b_values[23] * x7;
      d3 += b_values[31] * x7;
      d4 += b_values[39] * x7;
      d5 += b_values[47] * x7;
      d6 += b_values[55] * x7;
      d7 += b_values[63] * x7;
      d8 += b_values[71] * x7;
      d9 += b_values[79] * x7;
      d10 += b_values[87] * x7;
      d11 += b_values[95] * x7;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x9(const int &bm, const int *b_row_start, const int *b_col_idx,
               const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6, x7, x8;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 9) {
      x0 = x[9 * b_col_idx[j] + 0];
      x1 = x[9 * b_col_idx[j] + 1];
      x2 = x[9 * b_col_idx[j] + 2];
      x3 = x[9 * b_col_idx[j] + 3];
      x4 = x[9 * b_col_idx[j] + 4];
      x5 = x[9 * b_col_idx[j] + 5];
      x6 = x[9 * b_col_idx[j] + 6];
      x7 = x[9 * b_col_idx[j] + 7];
      x8 = x[9 * b_col_idx[j] + 8];
      d0 += b_values[0] * x0;
      d1 += b_values[9] * x0;
      d2 += b_values[18] * x0;
      d3 += b_values[27] * x0;
      d4 += b_values[36] * x0;
      d5 += b_values[45] * x0;
      d6 += b_values[54] * x0;
      d7 += b_values[63] * x0;
      d8 += b_values[72] * x0;
      d9 += b_values[81] * x0;
      d10 += b_values[90] * x0;
      d11 += b_values[99] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[10] * x1;
      d2 += b_values[19] * x1;
      d3 += b_values[28] * x1;
      d4 += b_values[37] * x1;
      d5 += b_values[46] * x1;
      d6 += b_values[55] * x1;
      d7 += b_values[64] * x1;
      d8 += b_values[73] * x1;
      d9 += b_values[82] * x1;
      d10 += b_values[91] * x1;
      d11 += b_values[100] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[11] * x2;
      d2 += b_values[20] * x2;
      d3 += b_values[29] * x2;
      d4 += b_values[38] * x2;
      d5 += b_values[47] * x2;
      d6 += b_values[56] * x2;
      d7 += b_values[65] * x2;
      d8 += b_values[74] * x2;
      d9 += b_values[83] * x2;
      d10 += b_values[92] * x2;
      d11 += b_values[101] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[12] * x3;
      d2 += b_values[21] * x3;
      d3 += b_values[30] * x3;
      d4 += b_values[39] * x3;
      d5 += b_values[48] * x3;
      d6 += b_values[57] * x3;
      d7 += b_values[66] * x3;
      d8 += b_values[75] * x3;
      d9 += b_values[84] * x3;
      d10 += b_values[93] * x3;
      d11 += b_values[102] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[13] * x4;
      d2 += b_values[22] * x4;
      d3 += b_values[31] * x4;
      d4 += b_values[40] * x4;
      d5 += b_values[49] * x4;
      d6 += b_values[58] * x4;
      d7 += b_values[67] * x4;
      d8 += b_values[76] * x4;
      d9 += b_values[85] * x4;
      d10 += b_values[94] * x4;
      d11 += b_values[103] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[14] * x5;
      d2 += b_values[23] * x5;
      d3 += b_values[32] * x5;
      d4 += b_values[41] * x5;
      d5 += b_values[50] * x5;
      d6 += b_values[59] * x5;
      d7 += b_values[68] * x5;
      d8 += b_values[77] * x5;
      d9 += b_values[86] * x5;
      d10 += b_values[95] * x5;
      d11 += b_values[104] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[15] * x6;
      d2 += b_values[24] * x6;
      d3 += b_values[33] * x6;
      d4 += b_values[42] * x6;
      d5 += b_values[51] * x6;
      d6 += b_values[60] * x6;
      d7 += b_values[69] * x6;
      d8 += b_values[78] * x6;
      d9 += b_values[87] * x6;
      d10 += b_values[96] * x6;
      d11 += b_values[105] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[16] * x7;
      d2 += b_values[25] * x7;
      d3 += b_values[34] * x7;
      d4 += b_values[43] * x7;
      d5 += b_values[52] * x7;
      d6 += b_values[61] * x7;
      d7 += b_values[70] * x7;
      d8 += b_values[79] * x7;
      d9 += b_values[88] * x7;
      d10 += b_values[97] * x7;
      d11 += b_values[106] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[17] * x8;
      d2 += b_values[26] * x8;
      d3 += b_values[35] * x8;
      d4 += b_values[44] * x8;
      d5 += b_values[53] * x8;
      d6 += b_values[62] * x8;
      d7 += b_values[71] * x8;
      d8 += b_values[80] * x8;
      d9 += b_values[89] * x8;
      d10 += b_values[98] * x8;
      d11 += b_values[107] * x8;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x10(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6, x7, x8, x9;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 10) {
      x0 = x[10 * b_col_idx[j] + 0];
      x1 = x[10 * b_col_idx[j] + 1];
      x2 = x[10 * b_col_idx[j] + 2];
      x3 = x[10 * b_col_idx[j] + 3];
      x4 = x[10 * b_col_idx[j] + 4];
      x5 = x[10 * b_col_idx[j] + 5];
      x6 = x[10 * b_col_idx[j] + 6];
      x7 = x[10 * b_col_idx[j] + 7];
      x8 = x[10 * b_col_idx[j] + 8];
      x9 = x[10 * b_col_idx[j] + 9];
      d0 += b_values[0] * x0;
      d1 += b_values[10] * x0;
      d2 += b_values[20] * x0;
      d3 += b_values[30] * x0;
      d4 += b_values[40] * x0;
      d5 += b_values[50] * x0;
      d6 += b_values[60] * x0;
      d7 += b_values[70] * x0;
      d8 += b_values[80] * x0;
      d9 += b_values[90] * x0;
      d10 += b_values[100] * x0;
      d11 += b_values[110] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[11] * x1;
      d2 += b_values[21] * x1;
      d3 += b_values[31] * x1;
      d4 += b_values[41] * x1;
      d5 += b_values[51] * x1;
      d6 += b_values[61] * x1;
      d7 += b_values[71] * x1;
      d8 += b_values[81] * x1;
      d9 += b_values[91] * x1;
      d10 += b_values[101] * x1;
      d11 += b_values[111] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[12] * x2;
      d2 += b_values[22] * x2;
      d3 += b_values[32] * x2;
      d4 += b_values[42] * x2;
      d5 += b_values[52] * x2;
      d6 += b_values[62] * x2;
      d7 += b_values[72] * x2;
      d8 += b_values[82] * x2;
      d9 += b_values[92] * x2;
      d10 += b_values[102] * x2;
      d11 += b_values[112] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[13] * x3;
      d2 += b_values[23] * x3;
      d3 += b_values[33] * x3;
      d4 += b_values[43] * x3;
      d5 += b_values[53] * x3;
      d6 += b_values[63] * x3;
      d7 += b_values[73] * x3;
      d8 += b_values[83] * x3;
      d9 += b_values[93] * x3;
      d10 += b_values[103] * x3;
      d11 += b_values[113] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[14] * x4;
      d2 += b_values[24] * x4;
      d3 += b_values[34] * x4;
      d4 += b_values[44] * x4;
      d5 += b_values[54] * x4;
      d6 += b_values[64] * x4;
      d7 += b_values[74] * x4;
      d8 += b_values[84] * x4;
      d9 += b_values[94] * x4;
      d10 += b_values[104] * x4;
      d11 += b_values[114] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[15] * x5;
      d2 += b_values[25] * x5;
      d3 += b_values[35] * x5;
      d4 += b_values[45] * x5;
      d5 += b_values[55] * x5;
      d6 += b_values[65] * x5;
      d7 += b_values[75] * x5;
      d8 += b_values[85] * x5;
      d9 += b_values[95] * x5;
      d10 += b_values[105] * x5;
      d11 += b_values[115] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[16] * x6;
      d2 += b_values[26] * x6;
      d3 += b_values[36] * x6;
      d4 += b_values[46] * x6;
      d5 += b_values[56] * x6;
      d6 += b_values[66] * x6;
      d7 += b_values[76] * x6;
      d8 += b_values[86] * x6;
      d9 += b_values[96] * x6;
      d10 += b_values[106] * x6;
      d11 += b_values[116] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[17] * x7;
      d2 += b_values[27] * x7;
      d3 += b_values[37] * x7;
      d4 += b_values[47] * x7;
      d5 += b_values[57] * x7;
      d6 += b_values[67] * x7;
      d7 += b_values[77] * x7;
      d8 += b_values[87] * x7;
      d9 += b_values[97] * x7;
      d10 += b_values[107] * x7;
      d11 += b_values[117] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[18] * x8;
      d2 += b_values[28] * x8;
      d3 += b_values[38] * x8;
      d4 += b_values[48] * x8;
      d5 += b_values[58] * x8;
      d6 += b_values[68] * x8;
      d7 += b_values[78] * x8;
      d8 += b_values[88] * x8;
      d9 += b_values[98] * x8;
      d10 += b_values[108] * x8;
      d11 += b_values[118] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[19] * x9;
      d2 += b_values[29] * x9;
      d3 += b_values[39] * x9;
      d4 += b_values[49] * x9;
      d5 += b_values[59] * x9;
      d6 += b_values[69] * x9;
      d7 += b_values[79] * x9;
      d8 += b_values[89] * x9;
      d9 += b_values[99] * x9;
      d10 += b_values[109] * x9;
      d11 += b_values[119] * x9;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x11(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6, x7, x8, x9, x10;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 11) {
      x0 = x[11 * b_col_idx[j] + 0];
      x1 = x[11 * b_col_idx[j] + 1];
      x2 = x[11 * b_col_idx[j] + 2];
      x3 = x[11 * b_col_idx[j] + 3];
      x4 = x[11 * b_col_idx[j] + 4];
      x5 = x[11 * b_col_idx[j] + 5];
      x6 = x[11 * b_col_idx[j] + 6];
      x7 = x[11 * b_col_idx[j] + 7];
      x8 = x[11 * b_col_idx[j] + 8];
      x9 = x[11 * b_col_idx[j] + 9];
      x10 = x[11 * b_col_idx[j] + 10];
      d0 += b_values[0] * x0;
      d1 += b_values[11] * x0;
      d2 += b_values[22] * x0;
      d3 += b_values[33] * x0;
      d4 += b_values[44] * x0;
      d5 += b_values[55] * x0;
      d6 += b_values[66] * x0;
      d7 += b_values[77] * x0;
      d8 += b_values[88] * x0;
      d9 += b_values[99] * x0;
      d10 += b_values[110] * x0;
      d11 += b_values[121] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[12] * x1;
      d2 += b_values[23] * x1;
      d3 += b_values[34] * x1;
      d4 += b_values[45] * x1;
      d5 += b_values[56] * x1;
      d6 += b_values[67] * x1;
      d7 += b_values[78] * x1;
      d8 += b_values[89] * x1;
      d9 += b_values[100] * x1;
      d10 += b_values[111] * x1;
      d11 += b_values[122] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[13] * x2;
      d2 += b_values[24] * x2;
      d3 += b_values[35] * x2;
      d4 += b_values[46] * x2;
      d5 += b_values[57] * x2;
      d6 += b_values[68] * x2;
      d7 += b_values[79] * x2;
      d8 += b_values[90] * x2;
      d9 += b_values[101] * x2;
      d10 += b_values[112] * x2;
      d11 += b_values[123] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[14] * x3;
      d2 += b_values[25] * x3;
      d3 += b_values[36] * x3;
      d4 += b_values[47] * x3;
      d5 += b_values[58] * x3;
      d6 += b_values[69] * x3;
      d7 += b_values[80] * x3;
      d8 += b_values[91] * x3;
      d9 += b_values[102] * x3;
      d10 += b_values[113] * x3;
      d11 += b_values[124] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[15] * x4;
      d2 += b_values[26] * x4;
      d3 += b_values[37] * x4;
      d4 += b_values[48] * x4;
      d5 += b_values[59] * x4;
      d6 += b_values[70] * x4;
      d7 += b_values[81] * x4;
      d8 += b_values[92] * x4;
      d9 += b_values[103] * x4;
      d10 += b_values[114] * x4;
      d11 += b_values[125] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[16] * x5;
      d2 += b_values[27] * x5;
      d3 += b_values[38] * x5;
      d4 += b_values[49] * x5;
      d5 += b_values[60] * x5;
      d6 += b_values[71] * x5;
      d7 += b_values[82] * x5;
      d8 += b_values[93] * x5;
      d9 += b_values[104] * x5;
      d10 += b_values[115] * x5;
      d11 += b_values[126] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[17] * x6;
      d2 += b_values[28] * x6;
      d3 += b_values[39] * x6;
      d4 += b_values[50] * x6;
      d5 += b_values[61] * x6;
      d6 += b_values[72] * x6;
      d7 += b_values[83] * x6;
      d8 += b_values[94] * x6;
      d9 += b_values[105] * x6;
      d10 += b_values[116] * x6;
      d11 += b_values[127] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[18] * x7;
      d2 += b_values[29] * x7;
      d3 += b_values[40] * x7;
      d4 += b_values[51] * x7;
      d5 += b_values[62] * x7;
      d6 += b_values[73] * x7;
      d7 += b_values[84] * x7;
      d8 += b_values[95] * x7;
      d9 += b_values[106] * x7;
      d10 += b_values[117] * x7;
      d11 += b_values[128] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[19] * x8;
      d2 += b_values[30] * x8;
      d3 += b_values[41] * x8;
      d4 += b_values[52] * x8;
      d5 += b_values[63] * x8;
      d6 += b_values[74] * x8;
      d7 += b_values[85] * x8;
      d8 += b_values[96] * x8;
      d9 += b_values[107] * x8;
      d10 += b_values[118] * x8;
      d11 += b_values[129] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[20] * x9;
      d2 += b_values[31] * x9;
      d3 += b_values[42] * x9;
      d4 += b_values[53] * x9;
      d5 += b_values[64] * x9;
      d6 += b_values[75] * x9;
      d7 += b_values[86] * x9;
      d8 += b_values[97] * x9;
      d9 += b_values[108] * x9;
      d10 += b_values[119] * x9;
      d11 += b_values[130] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[21] * x10;
      d2 += b_values[32] * x10;
      d3 += b_values[43] * x10;
      d4 += b_values[54] * x10;
      d5 += b_values[65] * x10;
      d6 += b_values[76] * x10;
      d7 += b_values[87] * x10;
      d8 += b_values[98] * x10;
      d9 += b_values[109] * x10;
      d10 += b_values[120] * x10;
      d11 += b_values[131] * x10;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}

void bcsr_12x12(const int &bm, const int *b_row_start, const int *b_col_idx,
                const double *b_values, const double *x, double *y) {
  int i, j;
  double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, x0, x1, x2, x3, x4,
      x5, x6, x7, x8, x9, x10, x11;
#pragma omp parallel for
  for (i = 0; i < bm; ++i) {
    d0 = y[12 * i + 0];
    d1 = y[12 * i + 1];
    d2 = y[12 * i + 2];
    d3 = y[12 * i + 3];
    d4 = y[12 * i + 4];
    d5 = y[12 * i + 5];
    d6 = y[12 * i + 6];
    d7 = y[12 * i + 7];
    d8 = y[12 * i + 8];
    d9 = y[12 * i + 9];
    d10 = y[12 * i + 10];
    d11 = y[12 * i + 11];
    for (j = b_row_start[i]; j < b_row_start[i + 1]; ++j, b_values += 12 * 12) {
      x0 = x[12 * b_col_idx[j] + 0];
      x1 = x[12 * b_col_idx[j] + 1];
      x2 = x[12 * b_col_idx[j] + 2];
      x3 = x[12 * b_col_idx[j] + 3];
      x4 = x[12 * b_col_idx[j] + 4];
      x5 = x[12 * b_col_idx[j] + 5];
      x6 = x[12 * b_col_idx[j] + 6];
      x7 = x[12 * b_col_idx[j] + 7];
      x8 = x[12 * b_col_idx[j] + 8];
      x9 = x[12 * b_col_idx[j] + 9];
      x10 = x[12 * b_col_idx[j] + 10];
      x11 = x[12 * b_col_idx[j] + 11];
      d0 += b_values[0] * x0;
      d1 += b_values[12] * x0;
      d2 += b_values[24] * x0;
      d3 += b_values[36] * x0;
      d4 += b_values[48] * x0;
      d5 += b_values[60] * x0;
      d6 += b_values[72] * x0;
      d7 += b_values[84] * x0;
      d8 += b_values[96] * x0;
      d9 += b_values[108] * x0;
      d10 += b_values[120] * x0;
      d11 += b_values[132] * x0;
      d0 += b_values[1] * x1;
      d1 += b_values[13] * x1;
      d2 += b_values[25] * x1;
      d3 += b_values[37] * x1;
      d4 += b_values[49] * x1;
      d5 += b_values[61] * x1;
      d6 += b_values[73] * x1;
      d7 += b_values[85] * x1;
      d8 += b_values[97] * x1;
      d9 += b_values[109] * x1;
      d10 += b_values[121] * x1;
      d11 += b_values[133] * x1;
      d0 += b_values[2] * x2;
      d1 += b_values[14] * x2;
      d2 += b_values[26] * x2;
      d3 += b_values[38] * x2;
      d4 += b_values[50] * x2;
      d5 += b_values[62] * x2;
      d6 += b_values[74] * x2;
      d7 += b_values[86] * x2;
      d8 += b_values[98] * x2;
      d9 += b_values[110] * x2;
      d10 += b_values[122] * x2;
      d11 += b_values[134] * x2;
      d0 += b_values[3] * x3;
      d1 += b_values[15] * x3;
      d2 += b_values[27] * x3;
      d3 += b_values[39] * x3;
      d4 += b_values[51] * x3;
      d5 += b_values[63] * x3;
      d6 += b_values[75] * x3;
      d7 += b_values[87] * x3;
      d8 += b_values[99] * x3;
      d9 += b_values[111] * x3;
      d10 += b_values[123] * x3;
      d11 += b_values[135] * x3;
      d0 += b_values[4] * x4;
      d1 += b_values[16] * x4;
      d2 += b_values[28] * x4;
      d3 += b_values[40] * x4;
      d4 += b_values[52] * x4;
      d5 += b_values[64] * x4;
      d6 += b_values[76] * x4;
      d7 += b_values[88] * x4;
      d8 += b_values[100] * x4;
      d9 += b_values[112] * x4;
      d10 += b_values[124] * x4;
      d11 += b_values[136] * x4;
      d0 += b_values[5] * x5;
      d1 += b_values[17] * x5;
      d2 += b_values[29] * x5;
      d3 += b_values[41] * x5;
      d4 += b_values[53] * x5;
      d5 += b_values[65] * x5;
      d6 += b_values[77] * x5;
      d7 += b_values[89] * x5;
      d8 += b_values[101] * x5;
      d9 += b_values[113] * x5;
      d10 += b_values[125] * x5;
      d11 += b_values[137] * x5;
      d0 += b_values[6] * x6;
      d1 += b_values[18] * x6;
      d2 += b_values[30] * x6;
      d3 += b_values[42] * x6;
      d4 += b_values[54] * x6;
      d5 += b_values[66] * x6;
      d6 += b_values[78] * x6;
      d7 += b_values[90] * x6;
      d8 += b_values[102] * x6;
      d9 += b_values[114] * x6;
      d10 += b_values[126] * x6;
      d11 += b_values[138] * x6;
      d0 += b_values[7] * x7;
      d1 += b_values[19] * x7;
      d2 += b_values[31] * x7;
      d3 += b_values[43] * x7;
      d4 += b_values[55] * x7;
      d5 += b_values[67] * x7;
      d6 += b_values[79] * x7;
      d7 += b_values[91] * x7;
      d8 += b_values[103] * x7;
      d9 += b_values[115] * x7;
      d10 += b_values[127] * x7;
      d11 += b_values[139] * x7;
      d0 += b_values[8] * x8;
      d1 += b_values[20] * x8;
      d2 += b_values[32] * x8;
      d3 += b_values[44] * x8;
      d4 += b_values[56] * x8;
      d5 += b_values[68] * x8;
      d6 += b_values[80] * x8;
      d7 += b_values[92] * x8;
      d8 += b_values[104] * x8;
      d9 += b_values[116] * x8;
      d10 += b_values[128] * x8;
      d11 += b_values[140] * x8;
      d0 += b_values[9] * x9;
      d1 += b_values[21] * x9;
      d2 += b_values[33] * x9;
      d3 += b_values[45] * x9;
      d4 += b_values[57] * x9;
      d5 += b_values[69] * x9;
      d6 += b_values[81] * x9;
      d7 += b_values[93] * x9;
      d8 += b_values[105] * x9;
      d9 += b_values[117] * x9;
      d10 += b_values[129] * x9;
      d11 += b_values[141] * x9;
      d0 += b_values[10] * x10;
      d1 += b_values[22] * x10;
      d2 += b_values[34] * x10;
      d3 += b_values[46] * x10;
      d4 += b_values[58] * x10;
      d5 += b_values[70] * x10;
      d6 += b_values[82] * x10;
      d7 += b_values[94] * x10;
      d8 += b_values[106] * x10;
      d9 += b_values[118] * x10;
      d10 += b_values[130] * x10;
      d11 += b_values[142] * x10;
      d0 += b_values[11] * x11;
      d1 += b_values[23] * x11;
      d2 += b_values[35] * x11;
      d3 += b_values[47] * x11;
      d4 += b_values[59] * x11;
      d5 += b_values[71] * x11;
      d6 += b_values[83] * x11;
      d7 += b_values[95] * x11;
      d8 += b_values[107] * x11;
      d9 += b_values[119] * x11;
      d10 += b_values[131] * x11;
      d11 += b_values[143] * x11;
      y[12 * i + 0] = d0;
      y[12 * i + 1] = d1;
      y[12 * i + 2] = d2;
      y[12 * i + 3] = d3;
      y[12 * i + 4] = d4;
      y[12 * i + 5] = d5;
      y[12 * i + 6] = d6;
      y[12 * i + 7] = d7;
      y[12 * i + 8] = d8;
      y[12 * i + 9] = d9;
      y[12 * i + 10] = d10;
      y[12 * i + 11] = d11;
    }
  }
}
