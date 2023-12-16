#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "headers/functions.hpp"
#include "headers/thomas.hpp"
#include "headers/de_allocate.hpp"

using namespace std;

/*
    Given a system of linear algebraic equations A * u = B,
    where A is a tridiagonal matrix such that
            β[0;N-1] is the main diagonal,
            α[0;N-2] is the diagonal over the main,
            γ[1;N-1] is the diagonal under the main;
          u is the column of unknowns (u[0;N-1]);
          B is the column of free terms (ω[0;N-1]).
    Set γ[0] and α[N-1] to zero.
*/

int main ()
{
    // Add test selection for single compilation
    string test;
    cout << "Choose a test!\nPossible input options: input3.txt, input4.txt, input5.txt" << '\n';
    cin >> test;
    cout << '\n';

    // Read from the file the number of equations and coefficients
    ifstream input(test);
    if (input.is_open())
    {
        int N = 0;
        input >> N;
        cout << "Hello, there're the solutions to the system of " << N << " equations for you!" << '\n';
    
        // Let's allocate memory for a two-dimensional array for an expanded matrix
        double ** A_expanded = create_array(N, N + 1);

        for (size_t i = 0; i != N; ++i)
            for (size_t j = 0; j != N + 1; ++j)
                input >> A_expanded[i][j];
        
        // Allocate memory for other necessary arrays
        double ** A = create_array(N, N);
        double * alpha = new double[N];
        double * beta = new double[N];
        double * gamma = new double[N];
        double * omega = new double[N];
        double * u = new double[N];

        // From the extended matrix, we get the matrix of coefficients, the column of free 
        // terms, and also fill in the arrays corresponding to the diagonals of the matrix
        for (size_t i = 0; i != N; ++i)
            for (size_t j = 0; j != N; ++j)
                A[i][j] = A_expanded[i][j];

        for (size_t i = 0; i != N; ++i)
        {
            beta[i] = A_expanded[i][i];
            omega[i] = A_expanded[i][N];
        }
                
        for (size_t i = 0; i != N - 1; ++i)
            alpha[i] = A_expanded[i][i + 1];
        alpha[N - 1] = 0;

        gamma[0] = 0;
        for (size_t i = 1; i != N; ++i)
            gamma[i] = A_expanded[i][i - 1];

        // Let's check Cramer's theorem
        const double eps = 1E-9;
        if (abs(trid_mat_determinant(N, alpha, beta, gamma)) < eps)
        {   // and check the Kronecker--Capelli theorem
            if (rk(N, N, A) != rk(N, N + 1, A_expanded))
                cout << "The system is unsolvable!" << '\n';
            else
                cout << "The solutions are endless!" << '\n';
        }
        else
        {   // Apply the tridiagonal matrix algorithm (Thomas algorithm)
            Thomas(N, alpha, beta, gamma, omega, u);
            
            // Printing the result to standard output
            for (size_t i = 0; i != N; ++i)
            {
                cout << round(u[i] * 1000) / 1000 << " ";
            }
            cout << '\n';

            // and also write the result to a file
            ofstream output("output.txt", ios::app);
            if (output.is_open())
                for (size_t i = 0; i != N; ++i)
                    output << round(u[i] * 1000) / 1000 << " ";
            output << '\n';
            output.close();  
        }
        
        // Finally, free the memory
        free_array(A_expanded);
        free_array(A);
        delete [] alpha;
        delete [] beta;
        delete [] gamma;
        delete [] omega;
        delete [] u;
    }
    input.close();

    return 0;
}
