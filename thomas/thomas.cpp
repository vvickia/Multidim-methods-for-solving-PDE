#include <iostream>

using namespace std;

// The function that implements the tridiagonal matrix algorithm
// (Derivation of formulas is presented in the report.)
void Thomas (int N, double *alpha, double *beta, double *gamma, double *omega, double *u)
{
    // Allocate the required memory
    double * x = new double[N];
    double * y = new double[N];

    // Reverse
    // by recursive formulas we find the values of x and y
    x[N - 2] = - gamma[N - 1] / beta[N - 1];
    y[N - 2] = omega[N - 1] / beta[N - 1];

    for (size_t j = N - 2; j >= 1; j--)
    {
        x[j - 1] = - gamma[j] / (alpha[j] * x[j] + beta[j]);
        y[j - 1] = (omega[j] - alpha[j] * y[j]) / (alpha[j] * x[j] + beta[j]);
    }

    // Direct move
    // thanks to the calculated initial value, we find the unknowns
    u[0] = (omega[0] - alpha[0] * y[0]) / (beta[0] + alpha[0] * x[0]);

    for (size_t j = 0; j != N - 1; ++j)
        u[j + 1] = x[j] * u[j] + y[j];

    // Freeing up memory
    delete [] x;
    delete [] y;
}
