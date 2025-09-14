#include "stdlib.h"
#include "string.h"
#include "matrix/clinalg_matrix.h"


void clinalg_matrix_set(clinalg_matrix_t *matrix, uint8_t row_index, uint8_t col_index, float value)
{
    matrix->data[row_index * matrix->cols + col_index] = value;
}

float clinalg_matrix_get(clinalg_matrix_t matrix, uint8_t row_index, uint8_t col_index)
{
    return matrix.data[row_index * matrix.cols + col_index];
}

void clinalg_matrix_set_data(clinalg_matrix_t *matrix, float *data)
{
    memcpy(matrix->data, data, sizeof(float) * matrix->rows * matrix->cols);
}

clinalg_matrix_t clinalg_matrix_add(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};

    // Assign matrix dimension
    matrix_result.rows = matrix_a.rows;
    matrix_result.cols = matrix_a.cols;

    for (uint8_t i = 0 ; i < matrix_a.rows ; i++)
    {
        for (uint8_t j = 0 ; j < matrix_a.cols ; j++)
        {
            matrix_result.data[i * matrix_a.cols + j] = matrix_a.data[i * matrix_a.cols + j] + matrix_b.data[i * matrix_a.cols + j];
        }
    }

    return matrix_result;
}


clinalg_matrix_t clinalg_matrix_subtract(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};

    // Assign matrix dimension
    matrix_result.rows = matrix_a.rows;
    matrix_result.cols = matrix_a.cols;

    for (uint8_t i = 0 ; i < matrix_a.rows ; i++)
    {
        for (uint8_t j = 0 ; j < matrix_a.cols ; j++)
        {
            matrix_result.data[i * matrix_a.cols + j] = matrix_a.data[i * matrix_a.cols + j] - matrix_b.data[i * matrix_a.cols + j];
        }
    }

    return matrix_result;
}


clinalg_matrix_t clinalg_matrix_multiply(clinalg_matrix_t matrix_a, clinalg_matrix_t matrix_b)
{
    clinalg_matrix_t matrix_result = {0};

    uint8_t m = matrix_a.rows;
    uint8_t n = matrix_a.cols;
    uint8_t k = matrix_b.cols;

    // Assign matrix dimension
    matrix_result.rows = m;
    matrix_result.cols = k;

    for (uint8_t i = 0 ; i < m ; i++)
    {
        for (uint8_t j = 0 ; j < k ; j++)
        {
            matrix_result.data[i * k + j] = matrix_a.data[i * n] * matrix_b.data[j];
            for (uint8_t s = 1; s < n ; s++)
            {
                matrix_result.data[i * k + j] += matrix_a.data[i * n + s] * matrix_b.data[s * k + j];
            }
        }
    }

    return matrix_result;
}

clinalg_matrix_t clinalg_matrix_transpose(clinalg_matrix_t matrix_a)
{
    clinalg_matrix_t matrix_result = {0};

    // Assign matrix dimension
    matrix_result.rows = matrix_a.cols;
    matrix_result.cols = matrix_a.rows;

    for (uint8_t i = 0 ; i < matrix_result.rows ; i++)
    {
        for (uint8_t j = 0 ; j < matrix_result.cols ; j++)
        {
            matrix_result.data[i * matrix_result.cols + j] = matrix_a.data[j * matrix_a.cols + i];
        }
    }

    return matrix_result;
}

float clinalg_matrix_determinant(clinalg_matrix_t matrix)
{
    uint8_t rows = matrix.rows;
    uint8_t cols = matrix.cols;

    if (rows == 1)
    {
        return matrix.data[0];
    }
    else if (rows == 2)
    {
        return matrix.data[0] * matrix.data[3] - matrix.data[1] * matrix.data[2];
    }

    float det = 0;
    int sign = 1;
    clinalg_matrix_t matrix_sub = {0};

    for (int col = 0; col < cols; col++)
    {
        int sub_index = 0;
        for (int i = 1; i < rows; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                if (j == col)
                    continue;
                matrix_sub.data[sub_index++] = matrix.data[i * rows + j];
            }
        }
        matrix_sub.rows = rows - 1;
        matrix_sub.cols = cols - 1;
        det += sign * matrix.data[col] * clinalg_matrix_determinant(matrix_sub);
        sign = -sign;
    }

    return det;
}

clinalg_matrix_t clinalg_matrix_cofactor(clinalg_matrix_t matrix, uint8_t p, uint8_t q)
{
    clinalg_matrix_t matrix_cofactor = {0};
    uint8_t N = matrix.rows;

    // Assign matrix dimension
    matrix_cofactor.rows = N - 1;
    matrix_cofactor.cols = N - 1;

    int i = 0, j = 0;
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            if (row != p && col != q)
            {
                matrix_cofactor.data[i * (N - 1) + j++] = matrix.data[row * N + col];
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
    clinalg_matrix_t matrix_cofactor = {0};

    uint8_t N = matrix.rows;
    if (N == 1)
    {
        matrix_adjugate.data[0] = 1;

        return matrix_adjugate;
    }

    int sign;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            matrix_cofactor = clinalg_matrix_cofactor(matrix, i, j);

            sign = ((i + j) % 2 == 0) ? 1 : -1;
            matrix_adjugate.data[j * N + i] = sign * clinalg_matrix_determinant(matrix_cofactor);
        }
    }

    return matrix_adjugate;
}

clinalg_matrix_t clinalg_matrix_inverse(clinalg_matrix_t matrix)
{
    clinalg_matrix_t matrix_inverse = {0};
    clinalg_matrix_t matrix_adjugate = {0};

    // Assign matrix dimension
    matrix_inverse.rows = matrix.rows;
    matrix_inverse.cols = matrix.cols;

    uint8_t N = matrix.rows;

    float det = clinalg_matrix_determinant(matrix);
    if (det == 0)
    {
        return matrix_inverse;
    }

    matrix_adjugate = clinalg_matrix_adjugate(matrix);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            matrix_inverse.data[i * N + j] = matrix_adjugate.data[i * N + j] / det;
        }
    }

    return matrix_inverse;
}


