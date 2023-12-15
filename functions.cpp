#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

// глобальные переменные (прописываю их и в hpp)
const double eps = 1e-8;
const int r = 100;
const double D = 1.11 * 1e-4;

void Initialization (double t_0, double t_b, double a, double b) // Task parameter initialization
{	
	ifstream inputFile("input.txt"); // Reading model parameters from a file
    if (inputFile.is_open())
    {
        inputFile >> a >> b >> t_0 >> t_b;
        inputFile.close();
    }
    else
    {
        cerr << "Ошибка открытия файла input.txt" << endl;
    }
}

void Thomas (vector<vector<double>> a, vector<double> b, vector<double> res)
{
    // real(8), dimension(0:, :) :: a
	// real(8), dimension(0:) :: b, res
	real(8), allocatable :: coef(:,:)
	real(8) eps
	integer i, lenght, ier
	
	lenght = size(b)
	
	// Memory allocation for array

	allocate(coef(0:lenght - 2, 0:1), stat = ier); if(ier /= 0) stop 4

	// X_N-2 and Y_N-2 calculation

	coef(lenght - 2, 0) = -a(lenght - 1, 1)/a(lenght - 1, 2)
	coef(lenght - 2, 1) = b(lenght - 1)/a(lenght - 1, 2)

	// Straight sweep

	do i = lenght - 3, 0, -1

		if (abs(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2)) < eps) then
		
			write(*,*) "Have not any results, or infinitely many results"; stop 10
			
		endif
		
		coef(i,0) = -a(i + 1, 1)/(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2))
		coef(i,1) = (b(i + 1) - (a(i + 1, 3) * coef(i + 1, 1)))/(a(i + 1, 3) * coef(i + 1, 0) + a(i + 1, 2))

	enddo
	
	// U_0 calculation
	if (abs(a(0,2) + a(0,3) * coef(0,0)) < eps) then
		
			write(*,*) "Have not any results, or infinitely many results"; stop 20
			
	endif
	res(0) = (b(0) - a(0,3) * coef(0,1))/ (a(0,2) + a(0,3) * coef(0,0))
	
	// Reverse sweep

	do i = 1, lenght - 1

		res(i) = coef(i - 1, 0) * res(i - 1) + coef(i - 1, 1)

	enddo

	deallocate(coef, stat = ier); if(ier /= 0) stop 5
}

void Set_IC (vector<vector<vector<double>>> u, double t_0, double t_b)
{
    for (size_t i = 1; i <= r + 2; ++i)
    {
        u[i][1][1] = t_b;
        u[i][1][2] = t_b;
        u[i][1][3] = t_b;
		u[i][r + 2][1] = t_b; // Border conditions
        u[i][r + 2][2] = t_b;
        u[i][r + 2][3] = t_b;
		
        for (size_t j = 2; j <= r + 1; ++j)
        {
            u[i][j][1] = t_0; // Initial
            u[i][j][2] = t_0;
            u[i][j][3] = t_0;
        }
    }
}

void Step (vector<vector<vector<double>>> u, double dt, double dx, double dy)
{
    vector<vector<double>> a;
    vector<double> b;
	
    a.resize(r, std::vector<double>(3));
    b.resize(r);
	
	// Stage 1

    for (size_t i = 1; i <= r - 2; ++i)
    {
        a[i][1] = -D / pow(dy, 2);
        a[i][3] = -D / pow(dy, 2);
		a[i][2] = 2 * D / pow(dy, 2) + 2 / dt;
    }
	
	a[0][1] = 0;
    a[0][2] = -1;
    a[0][3] = 1;
	a[r - 1][1] = -1;
    a[r - 1][2] = 1;
    a[r - 1][3] = 0;
	
	// a[0][1] = 0;
    // a[0][2] = 2 * D / pow(dy, 2) + 2 / dt;
    // a[0][3] = -D / pow(dy, 2);
	// a[r - 1][1] = -D / pow(dy, 2);
    // a[r - 1][2] = 2 * D / pow(dy, 2) + 2 / dt;
    // a[r - 1][3] = 0;
	
	for (size_t j = 2; j <= r + 1; ++j)
    {
        for (size_t i = 2; i <= r + 1; ++i)
        {
            b[i - 2] = D / pow(dx, 2) * (u[i + 1][j][1] - 2 * u[i][j][1] + u[i - 1][j][1]) + \
                       2 / dt * u[i][j][1];
        }
		
		b[0] = 0;
		b[r - 1] = 0;
		
		call Thomas(a, b, u(2:r+1,j,2)); // ???
    }
	
	Set_BC (u(:,:,2)); // Refresh BC
	
	
	// Stage 2
	
    for (size_t i = 1; i <= r - 2; ++i)
    {
        a[i][1] = -D / pow(dx, 2);
        a[i][3] = -D / pow(dx, 2);
		a[i][2] = 2 * D / pow(dx, 2) + 2 / dt;
    }

    a[0][1] = 0;
    a[0][2] = 1;
    a[0][3] = 0;
	a[r - 1][1] = 0;
    a[r - 1][2] = 1;
    a[r - 1][3] = 0;

    // a[0][1] = 0;
    // a[0][2] = 2 * D / pow(dx, 2) + 2 / dt;
    // a[0][3] = -D / pow(dx, 2);
	// a[r - 1][1] = -D / pow(dx, 2);
    // a[r - 1][2] = 2 * D / pow(dx, 2) + 2 / dt;
    // a[r - 1][3] = 0;

    for (size_t i = 2; i <= r + 1; ++i)
    {
        for (size_t j = 2; j <= r + 1; ++j)
        {
            b[j - 2] = D / pow(dy, 2) * (u[i][j + 1][2] - 2 * u[i][j][2] + u[i][j - 1][2]) + \
                       2 / dt * u[i][j][2];
        }
		
		b[0] = 25.0;
		b[r - 1] = 25.0;
		
		call Thomas(a, b, u(i,2:r+1,3)); // ???
    }
	
	Set_BC (u(:,:,3)); // Refresh BC
}

void Set_BC (vector<vector<double>> u)
{
    // u(1,:) = 0;
    // u(r+2,:) = 0
	
	u(1,:) = u(2,:);
    u(r + 2, :) = u(r + 1, :);
}

void Update_IC (vector<vector<vector<double>>> u)
{
    u(:,:,1) = u(:,:,3);
}

void Save_Data (vector<vector<vector<double>>> u, int iter)
{
    // Saving calculation results to a file.
	// Store the 2D matrix of results for each iteration step. 
	
    ostringstream str;
	
	if (iter < 10)
    {
        str << iter;
    }
    else if (iter < 100)
    {
        str << iter;
    }
    else
    {
        str << iter;
    }
	
    string filename = "result" + str.str() + ".txt";
    ofstream outputFile(filename);

    for (int i = 2; i <= r + 1; ++i)
    {
        write(1, *) u(i,:,3);
        // for (size_t j = 0; j < u[i][0].size(); ++j)
        // {
        //     outputFile << u[i][j][2] << " ";
        // }
        outputFile << "\n";
    }

    outputFile.close();
}
