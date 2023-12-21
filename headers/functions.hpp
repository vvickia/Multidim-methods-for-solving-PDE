#include <vector>

// using namespace std;

#ifndef PDE_FUNC
#define PDE_FUNC

extern const double eps = 1e-8;
extern const int r = 100;
extern const double D = 1.11 * 1e-4;

void Initialization(double &t_0, double &t_b, double &a, double &b);
void Thomas(const std::vector<std::vector<double>>& a, const std::vector<double>& b, std::vector<double>& res);
void Set_IC(std::vector<std::vector<std::vector<double>>>& u, double t_0, double t_b);
void Step(std::vector<std::vector<std::vector<double>>>& u, double dt, double dx, double dy);
void Set_BC(std::vector<std::vector<std::vector<double>>>& u, int layer);
void Update_IC(std::vector<std::vector<std::vector<double>>>& u);
void Save_Data(const std::vector<std::vector<std::vector<double>>>& u, int iter);

#endif