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

#define C_LINEAR_ALGEBRA_MAX_MATRIX_DATA 	32

typedef struct {
	uint8_t 	NumRows;
	uint8_t 	NumCols;
	float 		Data[C_LINEAR_ALGEBRA_MAX_MATRIX_DATA];
} CLinearAlgebra_Matrix_t;

void CLinearAlgebra_Matrix_SetValue(CLinearAlgebra_Matrix_t *Matrix, uint8_t RowIdx, uint8_t ColIdx, float Value);
float CLinearAlgebra_Matrix_GetValue(CLinearAlgebra_Matrix_t Matrix, uint8_t RowIdx, uint8_t ColIdx);
void CLinearAlgebra_Matrix_SetData(CLinearAlgebra_Matrix_t *Matrix, float *Data);
float CLinearAlgebra_Matrix_Determinant(CLinearAlgebra_Matrix_t Matrix);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Add(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Subtract(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Multiply(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Transpose(CLinearAlgebra_Matrix_t MatrixA);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Cofactor(CLinearAlgebra_Matrix_t Matrix, uint8_t p, uint8_t q);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Adjugate(CLinearAlgebra_Matrix_t Matrix);
CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Inverse(CLinearAlgebra_Matrix_t Matrix);
#ifdef __cplusplus
}
#endif

#endif /* __C_LINEAR_ALGEBRA_MATRIX_H__ */
