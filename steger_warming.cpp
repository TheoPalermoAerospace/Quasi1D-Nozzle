/*-----------------------------------------------------
------------- Theo Vinicius Souza Palermo -------------
-------------- Steger and Warming Scheme --------------
--------------- Quasi 1D Nozzle Problem ---------------
-----------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include "functions.h"

int main() {

	// p, p0, T, M and residual files
	std::ofstream dados_pressao("pressao.txt");
	std::ofstream dados_p0("pressao_estag.txt");
	std::ofstream dados_temp("temperatura.txt");
	std::ofstream dados_t0("temp_estag.txt");
	std::ofstream dados_mach("mach.txt");
	std::ofstream dados_r1("residuo1.txt");
	std::ofstream dados_r2("residuo2.txt");
	std::ofstream dados_r3("residuo3.txt");

	// Parametros
	const double gamma = 1.4;                      // Razao de Calores Especificos
	const double R = 287;                          // J/kg K
	const double CFL = 0.5;                        // Condicao CFL
	const int x_i = 0;                             // Coordenada do Primero Nó
	const int x_f = 1;                             // Coordenada do Ultimo Nó   
	const double L = 1.0;                          // Comprimento em x do Bocal CD
	const int N = 400;                             // Numero de Nos
	double delta_x = L/(N-1);                      // Espacamento da Malha
	double step = delta_x/2.0;
	std::vector<double> x = linspace(x_i, x_f, N);

	// Valores de referencia
	const double rho_ref = 1.486;                // kg/m3
	const double v_ref = 320.44;                 // m/s
	const double l_ref = 1.0;                    // m
	const double T_ref = 357.778;	              // K
	const double R_ref = 287;     	      // J/kg K

	// Constantes
	const double K1 = gamma/(gamma - 1);         // Termo do expoente de p_0
	const double K2 = (1 - 2*gamma)/(gamma - 1);
	
	// Condicoes iniciais de entrada
	double M_in = 0.7;      // Adimensional
	double p_in = 108.99e3; // Pa
	double T_in = 255.56;   // K

	// Variaveis primitivas e as demais: determinadas!
	double a_in = sqrt(gamma*R*T_in);                               // Velocidade do Som
	double u_in = M_in * a_in;                                      // Velocidade
	double rho_in = p_in/(R*T_in);                                  // Massa Especifica
	double e_in = (p_in / (gamma-1)) + 0.5*(rho_in * pow(u_in, 2)); // Energia Total

	// Condicoes iniciais de saida
	int caso;
	double M_out;
	double p_out;
	double T_out;
	double a_out;
	double e_out;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "             Caso 1     Caso 2 " << std::endl;
	std::cout << "-------------------------------" << std::endl;
	std::cout << "M_out [ - ]   0.70       0.70  " << std::endl;
	std::cout << "P_out [kPa]  85.000     40.000 " << std::endl;
	std::cout << "T_out [ K ]  255.56     255.56 " << std::endl;
	std::cout << "-------------------------------" << std::endl;
	std::cout << "\nInsira o caso desejado (1 ou 2): ";

	std::cin >> caso;

	if (caso == 1) {
		M_out = 0.7;                                               // Adimensional
		p_out = 85e3;                                              // Pa
		T_out = 255.56;                                            // K
		a_out = sqrt(gamma*R*T_out);                               // m/s
		e_out = (p_out / (gamma-1)) + 0.5*(rho_in * pow(u_in, 2)); // Energia Total
	}

	else if (caso == 2) {
		M_out = 0.7;                                               // Adimensional
		p_out = 40e3;                                              // Pa
		T_out = 255.56;                                            // K
		a_out = sqrt(gamma*R*T_out);                               // m/s
		e_out = (p_out / (gamma-1)) + 0.5*(rho_in * pow(u_in, 2)); // Energia Total
	}

	else {
		std::cout << "Entrada Invalida...\n" << std::endl;
	}

	// Inicializacao dos vetores
	std::vector<double> rho(N, 0.0);
	std::vector<double> u(N, 0.0);
	std::vector<double> a(N, 0.0);
	std::vector<double> a_carac(N, 0.0);
	std::vector<double> M(N, 0.0);
	std::vector<double> e(N, 0.0);
	std::vector<double> T(N, 0.0);
	std::vector<double> p(N, 0.0);
	std::vector<std::vector<double>> U(3, std::vector<double>(N, 0.0));   // Vetor de Estado
	std::vector<std::vector<double>> Q(3, std::vector<double>(N, 0.0));   // Termo Fonte
	std::vector<std::vector<double>> Res(3, std::vector<double>(N, 0.0)); // Residuos

	// Matrizes A_pos, A_neg, F_pos, F_neg, lambda_pos, lambda_neg, P e inv(P)
	std::vector<std::vector<double>> A_pos(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> A_neg(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> F_pos(3, std::vector<double> (N, 0.0));
	std::vector<std::vector<double>> F_neg(3, std::vector<double> (N, 0.0));
	std::vector<std::vector<double>> lambda_pos(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> lambda_neg(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> P(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> P_inv(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> U_aux(3, std::vector<double>(N, 0.0)); // Vetor auxiliar para a marcha no tempo
	
	// Vetores com os valores de area
	std::vector<double> S = calculate_area(x);
	std::vector<double> S_bar = calculate_sbar(x);

	// Adimensionalizacao
	double u_star = u_in/v_ref;
	double rho_star = rho_in/rho_ref;
	double p_in_star = p_in/(rho_ref * pow(v_ref, 2));
	double T_star = T_in/T_ref;
	double a_star = sqrt(gamma*T_star);
	double e_star = e_in/(rho_ref * pow(v_ref, 2));
	double e_star_out = e_out/(rho_ref * pow(v_ref, 2));
	double p_out_star = p_out/(rho_ref * pow(v_ref, 2));
	double R_star = R/R_ref;
	std::vector<double> x_star(x.size());

	for (size_t i = 0; i < x.size(); ++i) {
	    x_star[i] = x[i] / l_ref;
	}

	// Solucao Inicial: Escoamento Uniforme
	for (int i = 0; i < N; ++i) {
	    rho[i] = rho_star;
	    u[i] = u_star;
	    M[i] = M_in;
	    e[i] = e_star;
	    p[i] = p_in_star;
	    a[i] = a_star;
	    T[i] = T_star;
	}

	p[N-1] = p_out_star; // Valor de pressao do ultimo elemento
	e[N-1] = e_star_out; // Valor de energia total do ultimo elemento

	// Vetor de estado adimensionalizado
	for (int i = 0; i < N; ++i) {
		U[0][i] = rho[i];
		U[1][i] = rho[i] * u[i];
		U[2][i] = e[i];
	}

	U_aux = U;

	// Declaracao da pressao e da temperatura de estagnacao
	std::vector<double> p_0(N, 0.0);
	std::vector<double> T_0(N, 0.0);

	// Condicoes iniciais de p0 e T0
	for (int i = 0; i < N; ++i) {
		T_0[i] = T[i] * (1 + ( ((gamma - 1)/2)*pow(M[i], 2) ));
		p_0[i] = p[i] * pow((1 + ( ((gamma - 1)/2)*pow(M[i], 2) )), K1);
		a_carac[i] = sqrt(2*gamma*T_0[0] / (gamma+1));
	}

	// Declarando os termos R1, R2 e R3
	double R_1;
	double R_2;
	double R_3;

	// Declarando outras variaveis para a atualizacao
	double du_1;
	double dp_du;
	double du_N;
	double drho_N;
	double dp_N;

	// Declarando variáveis para R1, R2, R3
	double du_dx_1;
	double dp_dx_1;
	double dS_dx_1;
	double du_dx_N;
	double dp_dx_N;
	double dS_dx_N;
	double drho_dx_N;

	// Declarando os residuos
	double res_1;
	double res_2;
	double res_3;

	// Parametros da marcha temporal
	double criterio = 1e-6;  // Criterio de parada
	double residuo = 1;      // Apenas um valor inicial > 0
	int t_max = 20000;       // Instante máximo
	int t = 1;               // Instante inicial

	while (residuo > criterio && t < t_max) {

		// Calculo de delta_t
		double vel_max = maior_velocidade(U, a);
		double delta_t = CFL*delta_x/vel_max;

		/*-------------------------------------
		---------- Primeiro Elemento ----------
		-------------------------------------*/

		// Atualizar as variáveis no elemento i = 0
		du_dx_1 = (u[1] - u[0])/delta_x;
		dp_dx_1 = (p[1] - p[0])/delta_x;
		dS_dx_1 = (S[1] - S[0])/delta_x;

		R_3 = -(u[0] - a[0])*( du_dx_1 - (dp_dx_1/(rho[0]*a[0])) ) + ((u[0]*a[0]/S[0]) * dS_dx_1); // RHS da Equacao 2.19c 
		dp_du = -(p_0[0]*gamma*M[0]/a[0]) * pow((1 + 0.5*(gamma-1)*pow(M[0], 2)), K2);             // Equacao 5.13
		du_1 = (delta_t * R_3)/( 1 - (dp_du	/(rho[0]*a[0])) );                                     // Equacao 5.9
		
		u[0] += du_1;                                                                 // Equacao 5.14
		p[0] = p_0[0] * pow(1 - ((gamma-1)/(gamma+1))*pow(u[0]/a_carac[0], 2), K1);   // Equacao 5.15
		T[0] = T_0[0] * (1 - ((gamma-1)/(gamma+1))*pow(u[0]/a_carac[0], 2));          // Equacao 5.16
		rho[0] = p[0] / T[0];                                                         // Equacao 5.17
		e[0] = p[0]/(gamma-1) + (0.5*rho[0]*pow(u[0], 2));                            // Equacao 5.18
		a[0] = sqrt(gamma * T[0]);                                                    // Equacao 5.19
		M[0] = u[0]/a[0];                                                             // Equacao 5.20
		
		// Atualizacao do Vetor de Estado
		U[0][0] = rho[0];
		U[1][0] = rho[0] * u[0];
		U[2][0] = e[0];

		U_aux[0][0] = U[0][0];
		U_aux[1][0] = U[1][0];
		U_aux[2][0] = U[2][0];

		/*-----------------------------------
		---------- Ultimo Elemento ----------
		-----------------------------------*/
		du_dx_N = (u[N-1] - u[N-2])/delta_x;
		dp_dx_N = (p[N-1] - p[N-2])/delta_x;
		dS_dx_N = (S[N-1] - S[N-2])/delta_x;
		drho_dx_N = (rho[N-1] - rho[N-2])/delta_x;

		R_1 = -u[N-1] * ( drho_dx_N - (dp_dx_N/pow(a[N-1], 2)) );                                                // RHS da Equacao 2.19a
		R_2 = -(u[N-1] + a[N-2])*( du_dx_N + (dp_dx_N/(rho[N-1]*a[N-1])) ) - ((u[N-1]*a[N-1]/S[N-1]) * dS_dx_N); // RHS da Equacao 2.19b
		R_3 = -(u[N-1] - a[N-2])*( du_dx_N - (dp_dx_N/(rho[N-1]*a[N-1])) ) + ((u[N-1]*a[N-1]/S[N-1]) * dS_dx_N); // RHS da Equacao 2.19c

		// Avaliando numero de Mach no último elemento
		if (M[N-1] < 1) {

			p[N-1] = p_out_star;

			drho_N = R_1 * delta_t;                                      // Equacao 5.22
			du_N = R_2 * delta_t;                                        // Equacao 5.22

			rho[N-1] += drho_N;                                          // Equacao 5.23
			u[N-1] += du_N;                                              // Equacao 5.24
			e[N-1] = p[N-1]/(gamma-1) + ( 0.5*rho[N-1]*pow(u[N-1], 2) ); // Equacao 5.25
			T[N-1] = p[N-1]/rho[N-1];                                    // Equacao 5.26
			a[N-1] = sqrt(gamma * T[N-1]);                               // Equacao 5.27
			M[N-1] = u[N-1]/a[N-1];                                      // Equacao 5.28

			// Atualizacao do Vetor de Estado
			U[0][N-1] = rho[N-1];
			U[1][N-1] = rho[N-1] * u[N-1];
			U[2][N-1] = e[N-1];

			U_aux[0][N-1] = U[0][N-1];
			U_aux[1][N-1] = U[1][N-1];
			U_aux[2][N-1] = U[2][N-1];

		}

		else if (M[N-1] > 1) {

			du_N = (R_2 + R_3)*delta_t/2;                         // Equacao 5.29
			dp_N = rho[N-1] * a[N-1] * (R_2 - R_3) * delta_t / 2; // Equacao 5.30
			drho_N = (R_1 * delta_t) + (dp_N / pow(a[N-1], 2));   // Equacao 5.31

			u[N-1] += du_N;                                              // Equacao 5.32
			p[N-1] += dp_N;                                              // Equacao 5.33
			rho[N-1] += drho_N;                                          // Equacao 5.34
			e[N-1] = p[N-1]/(gamma-1) + ( 0.5*rho[N-1]*pow(u[N-1], 2) ); // Equacao 5.35
			T[N-1] = p[N-1]/rho[N-1];                                    // Equacao 5.36
			a[N-1] = sqrt(gamma * T[N-1]);                               // Equacao 5.37
			M[N-1] = u[N-1]/a[N-1];                                      // Equacao 5.38

			// Atualizacao do Vetor de Estado
			U[0][N-1] = rho[N-1];
			U[1][N-1] = rho[N-1] * u[N-1];
			U[2][N-1] = e[N-1];

			U_aux[0][N-1] = U[0][N-1];
			U_aux[1][N-1] = U[1][N-1];
			U_aux[2][N-1] = U[2][N-1];

		} else {

			std::cout << "Erro: numero de Mach negativo encontrado." << std::endl;

		}

		/*--------------------------------------
		---------- Elementos Internos ----------
		--------------------------------------*/

		// Nós internos (preenchimento das matrizes)
		for (int i = 0; i < N; ++i) {

			// Matriz de autovalores lambda_pos e lambda_neg
			if (M[i] < 1) {
				lambda_pos[0][0] = u[i];
				lambda_pos[1][1] = u[i] + a[i];
				lambda_pos[2][2] = 0.0;

				lambda_neg[0][0] = 0.0;
				lambda_neg[1][1] = 0.0;
				lambda_neg[2][2] = u[i] - a[i];
			
			} else if (M[i] > 1) {

				lambda_pos[0][0] = u[i];
				lambda_pos[1][1] = u[i] + a[i];
				lambda_pos[2][2] = u[i] - a[i];

				lambda_neg[0][0] = 0.0;
				lambda_neg[1][1] = 0.0;
				lambda_neg[2][2] = 0.0;
			}

			// Preencher P e inversa de P
			P[0][0] = 1.0;
			P[0][1] = U[0][i]/(2*a[i]);
			P[0][2] = -U[0][i]/(2*a[i]);
			P[1][0] = U[1][i]/U[0][i];
			P[1][1] = (U[0][i]/(2*a[i]))*((U[1][i]/U[0][i])+a[i]);
			P[1][2] = -(U[0][i]/(2*a[i]))*((U[1][i]/U[0][i])-a[i]);
			P[2][0] = pow((U[1][i]/U[0][i]), 2)/2;
			P[2][1] = (U[0][i]/(2*a[i]))*(0.5*pow((U[1][i]/U[0][i]), 2) + ((U[1][i]/U[0][i])*a[i]) + (pow(a[i], 2)/(gamma-1)));
			P[2][2] = -(U[0][i]/(2*a[i]))*(0.5*pow((U[1][i]/U[0][i]), 2) - ((U[1][i]/U[0][i])*a[i]) + (pow(a[i], 2)/(gamma-1)));

			P_inv[0][0] = 1 - (0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2)/pow(a[i], 2));
			P_inv[0][1] = (gamma-1)*(U[1][i]/U[0][i])/pow(a[i], 2);
			P_inv[0][2] = -(gamma-1)/pow(a[i], 2);
			P_inv[1][0] = (1/(U[0][i]*a[i]))*(0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2) - ((U[1][i]/U[0][i])*a[i]));
			P_inv[1][1] = (a[i] - ((gamma-1)*(U[1][i]/U[0][i])))/(U[0][i]*a[i]);
			P_inv[1][2] = (gamma-1)/(U[0][i]*a[i]);
			P_inv[2][0] = -(1/(U[0][i]*a[i]))*(0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2) + ((U[1][i]/U[0][i])*a[i]));
			P_inv[2][1] = (a[i] + ((gamma-1)*(U[1][i]/U[0][i])))/(U[0][i]*a[i]);
			P_inv[2][2] = -(gamma-1)/(U[0][i]*a[i]);

			// Jacobiano A_pos e A_neg
			std::vector<std::vector<double>> temp_pos = produto_matriz(P, lambda_pos); // P*lambda_pos
			A_pos = produto_matriz(temp_pos, P_inv);                                   // (P*lambda_pos)*P_inv
			
			std::vector<std::vector<double>> temp_neg = produto_matriz(P, lambda_neg); // P*lambda_neg
			A_neg = produto_matriz(temp_neg, P_inv);                                   // (P*lambda_neg)*P_inv

			// Preencher F_pos e F_neg
			std::vector<double> Fpos_i = produto_matriz_vetor(A_pos, {U[0][i], U[1][i], U[2][i]});        // F+_i
			std::vector<double> Fneg_i = produto_matriz_vetor(A_neg, {U[0][i], U[1][i], U[2][i]});        // F-_i
			std::vector<double> Fpos_im = produto_matriz_vetor(A_pos, {U[0][i-1], U[1][i-1], U[2][i-1]}); // F+_i+1
			std::vector<double> Fneg_ip = produto_matriz_vetor(A_neg, {U[0][i+1], U[1][i+1], U[2][i+1]}); // F-_i-1

			for (int j = 0; j < 3; ++j) {
	    		F_pos[j][i] = Fpos_i[j];
	    		F_neg[j][i] = Fneg_i[j];
			}

		}

		for (int i = 1; i < N-1; i++) {
			
			double S_ip = area((x[i] + x[i+1])/2.0); // S_{i+1/2}
			double S_im = area((x[i] + x[i-1])/2.0); // S_{i-1/2}
			//double Sbar_i = S_bar[i-1];              // S̄[i], pois S_bar vai de 0 até N-3
			double S_barra = (S_ip + S_im)/2;
			double dSdx = (S_ip - S_im)/delta_x;

			// Termo fonte Q_i^n
			Q[0][i] = 0.0;
			//Q[1][i] = (p[i]/S_bar[i-1]) * dS_dx[i-1];
			Q[1][i] = (p[i]/S_barra) * dSdx;
			Q[2][i] = 0.0;

            // Calculo dos residuos
			for (int j = 0; j < 3; ++j) {
				Res[j][i] = (-1.0 / S_barra) * ((S_ip*(F_pos[j][i] + F_neg[j][i+1]) - S_im*(F_pos[j][i-1] + F_neg[j][i])) / delta_x) + Q[j][i];
			}

			// Vetor de estado nos elementos internos
			for (int j = 0; j < 3; ++j) {
			    U_aux[j][i] = U[j][i] - (delta_t/(delta_x * S_barra)) * (S_ip*(F_pos[j][i]+F_neg[j][i+1]) - S_im*(F_pos[j][i-1]+F_neg[j][i])) + (delta_t*Q[j][i]);
			}

		}

		U = U_aux;

		for (int i = 0; i < N; ++i) {
            rho[i] = U[0][i];
            u[i] = U[1][i] / U[0][i];
            e[i] = U[2][i];
            p[i] = (e[i]-(0.5*rho[i]*pow(u[i], 2)))*(gamma-1);
            T[i] = p[i]/rho[i];
            a[i] = sqrt(gamma*T[i]);
            M[i] = (U[1][i]/U[0][i])/a[i];
            T_0[i] = T[i] * (1 + ((gamma - 1)/2)*pow(M[i], 2));
	        p_0[i] = p[i] * pow((1 + ( ((gamma - 1)/2)*pow(M[i], 2) )), K1);
        }

        std::cout << "Time Step: " << t << std::endl;

        residuo = calcular_residuos(Res, res_1, res_2, res_3);

        dados_r1 << t << " " << res_1 << std::endl;
        dados_r2 << t << " " << res_2 << std::endl;
        dados_r3 << t << " " << res_3 << std::endl;

		t++;

	}

	// Escrevendo nos arquivos
	for (int i = 0; i < N; ++i) {
    	dados_temp << x[i] << " " << T[i] << std::endl;
   		dados_p0 << x[i] << " " << p_0[i]  << std::endl;
   		dados_pressao << x[i] << " " << p[i]  << std::endl;
   		dados_mach << x[i] << " " << M[i]  << std::endl;
   		dados_t0 << x[i] << " " << T_0[i] << std::endl;
	}

    std::cout << "\n--------------------------------------" << std::endl;
    std::cout << "---------- Fim de Simulacao ----------" << std::endl;
    std::cout << "--------------------------------------" << std::endl;

    std::cout << "\nA simulação convergiu em " << t << " passos." << std::endl; 


    dados_pressao.close();
    dados_p0.close();
    dados_temp.close();
    dados_t0.close();
    dados_mach.close();
    dados_r1.close();
    dados_r2.close();
    dados_r3.close();
	
	return 0;
}
