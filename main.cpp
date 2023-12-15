#include <vector>

#include "headers/functions.hpp"

using namespace std;

int main ()
{
    double a, b;
    double t_0, t_b, dt, dx, dy;
    int n;

    vector<vector<vector<double>>> u;

    Initialization (t_0, t_b, a, b); // Task parameter initialization

    dx = a / r;
    dy = b / r;

    u.resize(r + 2, vector<vector<double>>(r + 2, vector<double>(3))); // создали 3-мерный вектор (r+2)x(r+2)x3

    Set_IC (u, t_0, t_b); // Setting initial conditions

    dt = 2; // Time step definition

    for (size_t i = 1; i <= 999; ++i)
    {
        Step (u, dt, dx, dy); // Calculation of the desired function at a new time step
        Update_IC (u); // Updating Initial Conditions
        Save_Data (u, i); // Saving calculation results to a file
    }

    return 0;
}

