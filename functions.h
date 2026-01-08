#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>

// Linspace Function
std::vector<double> linspace(double start, double end, int num);

// Nozzle Area
double area(double x);

// Calculate Area (S)
std::vector<double> calculate_area(const std::vector<double>& x);

// S_bar (\bar{S}) Calculation
std::vector<double> calculate_sbar(const std::vector<double>& x);

// Higher Characteristic Velocity
double maior_velocidade(std::vector<std::vector<double>>& U,
                        std::vector<double>& a);

// 3x3 Matrix Product
std::vector<std::vector<double>> produto_matriz(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B);

// Matrix-Vector Product
std::vector<double> produto_matriz_vetor(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& v);

// Vector (3x1) Summation
std::vector<double> soma_vetores(
    const std::vector<double>& a,
    const std::vector<double>& b);

// Residual Calculations
double calcular_residuos(const std::vector<std::vector<double>>& Res,
                         double& res_1,
                         double& res_2,
                         double& res_3);

#endif // FUNCTIONS_H
