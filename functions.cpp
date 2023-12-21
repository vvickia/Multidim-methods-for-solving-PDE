#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
// #include <sstream>

// using namespace std;

// глобальные переменные (прописываю их и в hpp)
const double eps = 1e-8;
const int r = 100;
const double D = 1.11 * 1e-4;

void Initialization(double &t_0, double &t_b, double &a, double &b)
{
    // Open file
    std::ifstream inputFile("input.txt");
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    // Read model parameters from file
    inputFile >> a >> b;
    inputFile >> t_0;
    inputFile >> t_b;

    // Close file
    inputFile.close();
}


void Thomas(const std::vector<std::vector<double>>& a, const std::vector<double>& b, std::vector<double>& res)
{
    int length = b.size();
    // double eps = 1e-10; // Define a small number for comparison with zero

    std::vector<std::vector<double>> coef(length, std::vector<double>(2));
    
    // X_N-2 and Y_N-2 calculation
    coef[length - 2][0] = -a[length - 1][1] / a[length - 1][2];
    coef[length - 2][1] = b[length - 1] / a[length - 1][2];

    // Straight sweep
    for (int i = length - 3; i >= 0; --i)
    {
        if (std::abs(a[i + 1][3] * coef[i + 1][0] + a[i + 1][2]) < eps)
        {
            std::cerr << "Have not any results, or infinitely many results" << std::endl;
            exit(10);
        }
        
        coef[i][0] = -a[i + 1][1] / (a[i + 1][3] * coef[i + 1][0] + a[i + 1][2]);
        coef[i][1] = (b[i + 1] - (a[i + 1][3] * coef[i + 1][1])) / (a[i + 1][3] * coef[i + 1][0] + a[i + 1][2]);
    }

    // U_0 calculation
    if (std::abs(a[0][2] + a[0][3] * coef[0][0]) < eps)
    {
        std::cerr << "Have not any results, or infinitely many results" << std::endl;
        exit(20);
    }
    res[0] = (b[0] - a[0][3] * coef[0][1]) / (a[0][2] + a[0][3] * coef[0][0]);

    // Reverse sweep
    for (int i = 1; i < length; ++i)
    {
        res[i] = coef[i - 1][0] * res[i - 1] + coef[i - 1][1];
    }
}


void Set_IC(std::vector<std::vector<std::vector<double>>>& u, double t_0, double t_b)
{
    // Assuming u is already initialized with the correct dimensions: (r + 2) x (r + 2) x 3

    for (int i = 0; i < r + 2; ++i)
    {
        // Set border conditions
        u[i][0][0] = t_b; u[i][0][1] = t_b; u[i][0][2] = t_b;
        u[i][r + 1][0] = t_b; u[i][r + 1][1] = t_b; u[i][r + 1][2] = t_b;

        // Set initial conditions for the inner values
        for (int j = 1; j <= r; ++j)
        {
            u[i][j][0] = t_0; u[i][j][1] = t_0; u[i][j][2] = t_0;
        }
    }
}

void Set_BC(std::vector<std::vector<std::vector<double>>>& u, int layer)
{
    // Обновление граничных условий для указанного слоя 'layer' в 3D векторе 'u'
    for (size_t j = 0; j < u[0].size(); ++j)
    {
        u[0][j][layer] = u[1][j][layer];
        u[r + 1][j][layer] = u[r][j][layer];
    }
}

void Step(std::vector<std::vector<std::vector<double>>>& u, double dt, double dx, double dy)
{
    int i, j;
    std::vector<std::vector<double>> a(r, std::vector<double>(3));
    std::vector<double> b(r - 1);

    // Stage 1
    for (i = 1; i < r - 1; ++i)
    {
        a[i][0] = -D / (dy * dy);
        a[i][2] = -D / (dy * dy);
        a[i][1] = 2 * D / (dy * dy) + 2 / dt;
    }

    a[0][0] = 0; a[0][1] = -1; a[0][2] = 1;
    a[r - 1][0] = -1; a[r - 1][1] = 1; a[r - 1][2] = 0;

    for (j = 1; j <= r; ++j)
    {
        for (i = 1; i <= r; ++i)
        {
            b[i - 1] = D / (dx * dx) * (u[i + 1][j][0] - 2 * u[i][j][0] + u[i - 1][j][0]) + 2 / dt * u[i][j][0];
        }
        b[0] = 0;
        b[r - 1] = 0;

        std::vector<double> res(r);
        Thomas(a, b, res);

        for (int k = 1; k <= r; ++k)
        {
            u[k][j][1] = res[k - 1];
        }
    }

    Set_BC(u, 1); // Refresh BC

    // Stage 2
    for (i = 1; i < r - 1; ++i)
    {
        a[i][0] = -D / (dx * dx);
        a[i][2] = -D / (dx * dx);
        a[i][1] = 2 * D / (dx * dx) + 2 / dt;
    }

    a[0][0] = 0; a[0][1] = 1; a[0][2] = 0;
    a[r - 1][0] = 0; a[r - 1][1] = 1; a[r - 1][2] = 0;

    for (i = 1; i <= r; ++i)
    {
        for (j = 1; j <= r; ++j)
        {
            b[j - 1] = D / (dy * dy) * (u[i][j + 1][1] - 2 * u[i][j][1] + u[i][j - 1][1]) + 2 / dt * u[i][j][1];
        }
        b[0] = 25.0; // Assuming these are the boundary conditions as per your code
        b[r - 1] = 25.0;

        std::vector<double> res(r);
        Thomas(a, b, res);

        for (int k = 1; k <= r; ++k)
        {
            u[i][k][2] = res[k - 1];
        }
    }

    Set_BC(u, 2); // Refresh BC

    // No dynamic memory to deallocate in C++
}

void Update_IC(std::vector<std::vector<std::vector<double>>>& u)
{
    // Assuming u is a 3D vector of dimensions [x][y][3]
    for (size_t i = 0; i < u.size(); ++i)
    {
        for (size_t j = 0; j < u[i].size(); ++j)
        {
            u[i][j][0] = u[i][j][2]; // Copying the third layer to the first layer
        }
    }
}

void Save_Data(const std::vector<std::vector<std::vector<double>>>& u, int iter)
{
    // Format the file name based on the iteration number
    std::string filename = "result" + std::to_string(iter) + ".dat";

    // Create and open a file stream
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    // Write the 2D matrix (third layer of u) to the file
    for (int i = 1; i <= r; ++i)
    {
        for (int j = 0; j < u[i].size(); ++j)
        {
            file << std::setprecision(8) << std::fixed << u[i][j][2];
            if (j < u[i].size() - 1)
            {
                file << " "; // Add space between numbers in the same row
            }
        }
        file << "\n"; // New line for each row
    }

    // Close the file
    file.close();
}

