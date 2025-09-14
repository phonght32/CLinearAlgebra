// MIT License

// Copyright (c) 2025 phonght32

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef __C_LINEAR_ALGEBRA_MATRIX_H__
#define __C_LINEAR_ALGEBRA_MATRIX_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "stdint.h"

#ifndef C_LINEAR_ALGEBRA_MAX_ROWS
#define C_LINEAR_ALGEBRA_MAX_ROWS   4
#endif

#ifndef C_LINEAR_ALGEBRA_MAX_COLS
#define C_LINEAR_ALGEBRA_MAX_COLS   4
#endif

typedef struct
{
    uint8_t     rows;
    uint8_t     cols;
    float       data[C_LINEAR_ALGEBRA_MAX_ROWS][C_LINEAR_ALGEBRA_MAX_COLS];
} clinalg_matrix_t;

void clinalg_matrix_set(clinalg_matrix_t *matrix, uint8_t row_index, uint8_t col_index, float Value);
float clinalg_matrix_get(clinalg_matrix_t matrix, uint8_t row_index, uint8_t col_index);
void clinalg_matrix_set_data(clinalg_matrix_t *matrix, float *data);
float clinalg_matrix_determinant(clinalg_matrix_t matrix);
clinalg_matrix_t clinalg_matrix_add(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b);
clinalg_matrix_t clinalg_matrix_subtract(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b);
clinalg_matrix_t clinalg_matrix_multiply(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b);
clinalg_matrix_t clinalg_matrix_transpose(clinalg_matrix_t matrix_a);
clinalg_matrix_t clinalg_matrix_cofactor(clinalg_matrix_t matrix, uint8_t p, uint8_t q);
clinalg_matrix_t clinalg_matrix_adjugate(clinalg_matrix_t matrix);
clinalg_matrix_t clinalg_matrix_inverse(clinalg_matrix_t matrix);
#ifdef __cplusplus
}
#endif

#endif /* __C_LINEAR_ALGEBRA_MATRIX_H__ */
