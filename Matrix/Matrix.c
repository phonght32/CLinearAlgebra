#include "stdlib.h"
#include "string.h"

#include "Matrix/Matrix.h"

//static void SwapRows(CLinearAlgebra_Matrix_t *Matrix, uint8_t RowIdx1, uint8_t RowIdx2)
//{
//    if ((Matrix->NumRows <= RowIdx1) || (Matrix->NumRows <= RowIdx2))
//    {
//        return;
//    }
//    else
//    {
//        uint8_t stride = Matrix->NumCols;
//        for (int i = 0; i < Matrix->NumCols; i++)
//        {
//            float temp = Matrix->Data[RowIdx1 * stride + i];
//            Matrix->Data[RowIdx1 * stride + i] = Matrix->Data[RowIdx2 * stride + i];
//            Matrix->Data[RowIdx2 * stride + i] = temp;
//        }
//    }
//}

void CLinearAlgebra_Matrix_SetValue(CLinearAlgebra_Matrix_t *Matrix, uint8_t RowIdx, uint8_t ColIdx, float Val)
{
    Matrix->Data[RowIdx * Matrix->NumCols + ColIdx] = Val;
}

float CLinearAlgebra_Matrix_GetValue(CLinearAlgebra_Matrix_t Matrix, uint8_t RowIdx, uint8_t ColIdx)
{
    return Matrix.Data[RowIdx * Matrix.NumCols + ColIdx];
}

void CLinearAlgebra_Matrix_SetData(CLinearAlgebra_Matrix_t *Matrix, float *Data)
{
    memcpy(Matrix->Data, Data, sizeof(float) * Matrix->NumRows * Matrix->NumCols);
}

CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Add(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB)
{
    CLinearAlgebra_Matrix_t MatrixResult = {0};

    // Assign matrix dimension
    MatrixResult.NumRows = MatrixA.NumRows;
    MatrixResult.NumCols = MatrixA.NumCols;

    for (uint8_t i = 0 ; i < MatrixA.NumRows ; i++)
    {
        for (uint8_t j = 0 ; j < MatrixA.NumCols ; j++)
        {
            MatrixResult.Data[i * MatrixA.NumCols + j] = MatrixA.Data[i * MatrixA.NumCols + j] + MatrixB.Data[i * MatrixA.NumCols + j];
        }
    }

    return MatrixResult;
}


CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Subtract(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB)
{
    CLinearAlgebra_Matrix_t MatrixResult = {0};

    // Assign matrix dimension
    MatrixResult.NumRows = MatrixA.NumRows;
    MatrixResult.NumCols = MatrixA.NumCols;

    for (uint8_t i = 0 ; i < MatrixA.NumRows ; i++)
    {
        for (uint8_t j = 0 ; j < MatrixA.NumCols ; j++)
        {
            MatrixResult.Data[i * MatrixA.NumCols + j] = MatrixA.Data[i * MatrixA.NumCols + j] - MatrixB.Data[i * MatrixA.NumCols + j];
        }
    }

    return MatrixResult;
}


CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Multiply(CLinearAlgebra_Matrix_t MatrixA, CLinearAlgebra_Matrix_t MatrixB)
{
    CLinearAlgebra_Matrix_t MatrixResult = {0};

    uint8_t m = MatrixA.NumRows;
    uint8_t n = MatrixA.NumCols;
    uint8_t k = MatrixB.NumCols;

    // Assign matrix dimension
    MatrixResult.NumRows = m;
    MatrixResult.NumCols = k;

    for (uint8_t i = 0 ; i < m ; i++)
    {
        for (uint8_t j = 0 ; j < k ; j++)
        {
            MatrixResult.Data[i * k + j] = MatrixA.Data[i * n] * MatrixB.Data[j];
            for (uint8_t s = 1; s < n ; s++)
            {
                MatrixResult.Data[i * k + j] += MatrixA.Data[i * n + s] * MatrixB.Data[s * k + j];
            }
        }
    }

    return MatrixResult;
}

CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Transpose(CLinearAlgebra_Matrix_t MatrixA)
{
    CLinearAlgebra_Matrix_t MatrixResult = {0};

    // Assign matrix dimension
    MatrixResult.NumRows = MatrixA.NumCols;
    MatrixResult.NumCols = MatrixA.NumRows;

    for (uint8_t i = 0 ; i < MatrixResult.NumRows ; i++)
    {
        for (uint8_t j = 0 ; j < MatrixResult.NumCols ; j++)
        {
            MatrixResult.Data[i * MatrixResult.NumCols + j] = MatrixA.Data[j * MatrixA.NumCols + i];
        }
    }

    return MatrixResult;
}

float CLinearAlgebra_Matrix_Determinant(CLinearAlgebra_Matrix_t Matrix)
{
    uint8_t NumRows = Matrix.NumRows;
    uint8_t NumCols = Matrix.NumCols;

    if (NumRows == 1)
    {
        return Matrix.Data[0];
    }
    else if (NumRows == 2)
    {
        return Matrix.Data[0] * Matrix.Data[3] - Matrix.Data[1] * Matrix.Data[2];
    }

    float det = 0;
    int sign = 1;
    CLinearAlgebra_Matrix_t SubMatrix = {0};

    for (int col = 0; col < NumCols; col++)
    {
        int subIndex = 0;
        for (int i = 1; i < NumRows; i++)
        {
            for (int j = 0; j < NumRows; j++)
            {
                if (j == col)
                    continue;
                SubMatrix.Data[subIndex++] = Matrix.Data[i * NumRows + j];
            }
        }
        SubMatrix.NumRows = NumRows - 1;
        SubMatrix.NumCols = NumCols - 1;
        det += sign * Matrix.Data[col] * CLinearAlgebra_Matrix_Determinant(SubMatrix);
        sign = -sign;
    }

    return det;
}

CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Cofactor(CLinearAlgebra_Matrix_t Matrix, uint8_t p, uint8_t q)
{
    CLinearAlgebra_Matrix_t MatrixCofactor = {0};
    uint8_t N = Matrix.NumRows;

    // Assign matrix dimension
    MatrixCofactor.NumRows = N - 1;
    MatrixCofactor.NumCols = N - 1;

    int i = 0, j = 0;
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            if (row != p && col != q)
            {
                MatrixCofactor.Data[i * (N - 1) + j++] = Matrix.Data[row * N + col];
                if (j == N - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return MatrixCofactor;
}

CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Adjugate(CLinearAlgebra_Matrix_t Matrix)
{
    CLinearAlgebra_Matrix_t MatrixAdjugate = {0};
    CLinearAlgebra_Matrix_t MatrixCofactor = {0};

    uint8_t N = Matrix.NumRows;
    if (N == 1)
    {
        MatrixAdjugate.Data[0] = 1;

        return MatrixAdjugate;
    }

    int sign;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            MatrixCofactor = CLinearAlgebra_Matrix_Cofactor(Matrix, i, j);

            sign = ((i + j) % 2 == 0) ? 1 : -1;
            MatrixAdjugate.Data[j * N + i] = sign * CLinearAlgebra_Matrix_Determinant(MatrixCofactor);
        }
    }

    return MatrixAdjugate;
}

CLinearAlgebra_Matrix_t CLinearAlgebra_Matrix_Inverse(CLinearAlgebra_Matrix_t Matrix)
{
    CLinearAlgebra_Matrix_t MatrixInverse = {0};
    CLinearAlgebra_Matrix_t MatrixAdjugate = {0};

    // Assign matrix dimension
    MatrixInverse.NumRows = Matrix.NumRows;
    MatrixInverse.NumCols = Matrix.NumCols;

    uint8_t N = Matrix.NumRows;

    float det = CLinearAlgebra_Matrix_Determinant(Matrix);
    if (det == 0)
    {
        return MatrixInverse;
    }

    MatrixAdjugate = CLinearAlgebra_Matrix_Adjugate(Matrix);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            MatrixInverse.Data[i * N + j] = MatrixAdjugate.Data[i * N + j] / det;
        }
    }

    return MatrixInverse;
}


