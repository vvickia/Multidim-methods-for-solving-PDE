#include <vector>

using namespace std;

#ifndef PDE_FUNC
#define PDE_FUNC

extern const double eps = 1e-8;
extern const int r = 100;
extern const double D = 1.11 * 1e-4;

void Initialization (double t_0, double t_b, double a, double b);
void Thomas (vector<vector<double>> a, vector<double> b, vector<double> res);
void Set_IC (vector<vector<vector<double>>> u, double t_0, double t_b);
void Step (vector<vector<vector<double>>>u, double dt, double dx, double dy);
void Set_BC (vector<vector<double>> u);
void Update_IC (vector<vector<vector<double>>> u);
void Save_Data (vector<vector<vector<double>>> u, int iter);

#endif