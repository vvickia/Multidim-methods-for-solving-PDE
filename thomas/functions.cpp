#include <cstdlib>
#include <vector>

using namespace std;

// The function of finding the determinant of a tridiagonal matrix
double trid_mat_determinant (int N, double *alpha, double *beta, double *gamma)
{
    double det_prev_prev = abs(beta[0]);
    double det_prev = beta[1] * det_prev_prev - gamma[1] * alpha[0];
    double det = 0; // doesn't work with 2x2 matrices
    double tmp = 0;

    for (size_t i = 2; i != N; ++i)
    {
        det = beta[i] * det_prev - gamma[i] * alpha[i - 1] * det_prev_prev;
        tmp = det_prev;
        det_prev = det;
        det_prev_prev = tmp;
    }

    return det;
}

// The function of finding the rank of a matrix
int rk (int rows, int cols, double **mat)
{
    const double eps = 1E-9;
    int rank = cols;

    vector<char> line_used (rows);
    for (int i = 0; i < cols; ++i)
    {
        int j;
        for (j = 0; j < rows; ++j)
            if (!line_used[j] && abs(mat[j][i]) > eps)
                break;
        if (j == rows)
            --rank;
        else
        {
            line_used[j] = true;
            for (int p = i + 1; p < cols; ++p)
                mat[j][p] /= mat[j][i];
            for (int k = 0; k < rows; ++k)
                if (k != j && abs(mat[k][i]) > eps)
                    for (int p = i + 1; p < cols; ++p)
                        mat[k][p] -= mat[j][p] * mat[k][i];
        }
    }
    return rank;
}
