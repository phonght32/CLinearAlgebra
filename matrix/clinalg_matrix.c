#include "stdlib.h"
#include "string.h"
#include "matrix/clinalg_matrix.h"

clinalg_matrix_t clinalg_matrix_create(uint8_t rows, uint8_t cols, float *data)
{
    clinalg_matrix_t matrix = {0};
    matrix.rows = rows;
    matrix.cols = cols;

    for (uint8_t i = 0; i < rows; i++)
    {
        for (uint8_t j = 0; j < cols; j++)
        {
            matrix.data[i][j] = data[i * cols + j];
        }
    }

    return matrix;
}

void clinalg_matrix_set_member(clinalg_matrix_t *matrix, uint8_t row_index, uint8_t col_index, float value)
{
    matrix->data[row_index][col_index] = value;
}

float clinalg_matrix_get_member(clinalg_matrix_t matrix, uint8_t row_index, uint8_t col_index)
{
    return matrix.data[row_index][col_index];
}

clinalg_matrix_t clinalg_matrix_add(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};
    matrix_result.rows = matrix_a.rows;
    matrix_result.cols = matrix_a.cols;

    for (uint8_t i = 0; i < matrix_a.rows; i++)
    {
        for (uint8_t j = 0; j < matrix_a.cols; j++)
        {
            matrix_result.data[i][j] = matrix_a.data[i][j] + matrix_b.data[i][j];
        }
    }

    return matrix_result;
}

clinalg_matrix_t clinalg_matrix_subtract(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};
    matrix_result.rows = matrix_a.rows;
    matrix_result.cols = matrix_a.cols;

    for (uint8_t i = 0; i < matrix_a.rows; i++)
    {
        for (uint8_t j = 0; j < matrix_a.cols; j++)
        {
            matrix_result.data[i][j] = matrix_a.data[i][j] - matrix_b.data[i][j];
        }
    }

    return matrix_result;
}

clinalg_matrix_t clinalg_matrix_multiply(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};
    matrix_result.rows = matrix_a.rows;
    matrix_result.cols = matrix_b.cols;

    for (uint8_t i = 0; i < matrix_a.rows; i++)
    {
        for (uint8_t j = 0; j < matrix_b.cols; j++)
        {
            matrix_result.data[i][j] = 0;
            for (uint8_t k = 0; k < matrix_a.cols; k++)
            {
                matrix_result.data[i][j] += matrix_a.data[i][k] * matrix_b.data[k][j];
            }
        }
    }

    return matrix_result;
}

clinalg_matrix_t clinalg_matrix_transpose(clinalg_matrix_t matrix_a)
{
    clinalg_matrix_t matrix_result = {0};
    matrix_result.rows = matrix_a.cols;
    matrix_result.cols = matrix_a.rows;

    for (uint8_t i = 0; i < matrix_result.rows; i++)
    {
        for (uint8_t j = 0; j < matrix_result.cols; j++)
        {
            matrix_result.data[i][j] = matrix_a.data[j][i];
        }
    }

    return matrix_result;
}

float clinalg_matrix_determinant(clinalg_matrix_t matrix)
{
    if (matrix.rows == 1)
        return matrix.data[0][0];
    if (matrix.rows == 2)
        return matrix.data[0][0] * matrix.data[1][1] - matrix.data[0][1] * matrix.data[1][0];

    float det = 0;
    int sign = 1;
    clinalg_matrix_t sub = {0};
    sub.rows = matrix.rows - 1;
    sub.cols = matrix.cols - 1;

    for (int col = 0; col < matrix.cols; col++)
    {
        int subi = 0;
        for (int i = 1; i < matrix.rows; i++)
        {
            int subj = 0;
            for (int j = 0; j < matrix.cols; j++)
            {
                if (j == col) continue;
                sub.data[subi][subj++] = matrix.data[i][j];
            }
            subi++;
        }

        det += sign * matrix.data[0][col] * clinalg_matrix_determinant(sub);
        sign = -sign;
    }

    return det;
}

clinalg_matrix_t clinalg_matrix_cofactor(clinalg_matrix_t matrix, uint8_t p, uint8_t q)
{
    clinalg_matrix_t matrix_cofactor = {0};
    uint8_t N = matrix.rows;
    matrix_cofactor.rows = N - 1;
    matrix_cofactor.cols = N - 1;

    int i = 0, j = 0;
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            if (row != p && col != q)
            {
                matrix_cofactor.data[i][j++] = matrix.data[row][col];
                if (j == N - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return matrix_cofactor;
}

clinalg_matrix_t clinalg_matrix_adjugate(clinalg_matrix_t matrix)
{
    clinalg_matrix_t matrix_adjugate = {0};
    uint8_t N = matrix.rows;
    matrix_adjugate.rows = N;
    matrix_adjugate.cols = N;

    if (N == 1)
    {
        matrix_adjugate.data[0][0] = 1;
        return matrix_adjugate;
    }

    int sign;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            clinalg_matrix_t cof = clinalg_matrix_cofactor(matrix, i, j);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            matrix_adjugate.data[j][i] = sign * clinalg_matrix_determinant(cof);
        }
    }

    return matrix_adjugate;
}

clinalg_matrix_t clinalg_matrix_inverse(clinalg_matrix_t matrix)
{
    clinalg_matrix_t matrix_inverse = {0};
    matrix_inverse.rows = matrix.rows;
    matrix_inverse.cols = matrix.cols;

    float det = clinalg_matrix_determinant(matrix);
    if (det == 0)
        return matrix_inverse;

    clinalg_matrix_t adj = clinalg_matrix_adjugate(matrix);

    for (int i = 0; i < matrix.rows; i++)
    {
        for (int j = 0; j < matrix.cols; j++)
        {
            matrix_inverse.data[i][j] = adj.data[i][j] / det;
        }
    }

    return matrix_inverse;
}
