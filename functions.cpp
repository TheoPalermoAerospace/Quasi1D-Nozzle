#include "functions.h"

#include <cmath>
#include <algorithm>

// Linspace Functions
std::vector<double> linspace(double start, double end, int num) {

    std::vector<double> x;
    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        x.push_back(start + i * step);
    }

    return x;
}

// Nozzle Area
double area(double x) {
    const double a = 0.8;
    const double b = -0.6;
    const double c = 0.35;

    return a * x * x + b * x + c;
}

// Calculate Area (S)
std::vector<double> calculate_area(const std::vector<double>& x) {

    std::vector<double> S;

    for (size_t i = 0; i < x.size(); ++i) {
        S.push_back(area(x[i]));
    }

    return S;
}

// S_bar (\bar{S}) Calculation
std::vector<double> calculate_sbar(const std::vector<double>& x) {

    std::vector<double> S_bar;

    for (size_t i = 1; i < x.size() - 1; ++i) {
        double x_i_menos = (x[i] + x[i - 1]) / 2.0;
        double x_i_mais  = (x[i] + x[i + 1]) / 2.0;

        S_bar.push_back((area(x_i_menos) + area(x_i_mais)) / 2.0);
    }

    return S_bar;
}

// Higher Characteristic Velocity
double maior_velocidade(std::vector<std::vector<double>>& U,
                        std::vector<double>& a) {

    double vel_max = 0.0;
    size_t n = U[0].size();

    for (size_t i = 0; i < n; ++i) {

        double densidade = U[0][i];
        double velocidade = U[1][i] / densidade;

        double valor = std::abs(velocidade + a[i]);
        vel_max = std::max(vel_max, valor);
    }

    return vel_max;
}

// 3x3 Matrix Product
std::vector<std::vector<double>> produto_matriz(
            const std::vector<std::vector<double>>& A,
            const std::vector<std::vector<double>>& B) {

    std::vector<std::vector<double>> produto(3, std::vector<double>(3, 0.0));

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                produto[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return produto;
}

// Matrix-Vector Product
std::vector<double> produto_matriz_vetor(
            const std::vector<std::vector<double>>& A,
            const std::vector<double>& v) {

    std::vector<double> prod(3, 0.0);

    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
            prod[i] += A[i][k] * v[k];
        }
    }

    return prod;
}

// Vector (3x1) Summation
std::vector<double> soma_vetores(
            const std::vector<double>& a,
            const std::vector<double>& b) {

    std::vector<double> res(3);

    for (int k = 0; k < 3; ++k) {
        res[k] = a[k] + b[k];
    }

    return res;
}

// Residual Calculations
double calcular_residuos(const std::vector<std::vector<double>>& Res,
                         double& res_1,
                         double& res_2,
                         double& res_3) {

    res_1 = 0.0;
    res_2 = 0.0;
    res_3 = 0.0;

    int N = Res[0].size();

    for (int i = 1; i < N - 1; ++i) {
        res_1 = std::max(res_1, std::abs(Res[0][i]));
        res_2 = std::max(res_2, std::abs(Res[1][i]));
        res_3 = std::max(res_3, std::abs(Res[2][i]));
    }

    return std::max({res_1, res_2, res_3});
}
