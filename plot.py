#-------------------------------------------------------#
#------------- Theo Vinicius Souza Palermo -------------#
#----- Dinamica dos Fluidos Computacional: Q2.2025 -----#
#-------------- Trabalho Final - Bocal CD --------------#
#-------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt

# Funcao para a leitura dos arquivos
def ler_dados(nome_arquivo):
    return np.loadtxt(nome_arquivo)

# Leitura dos arquivos
x_p, p = ler_dados("pressao.txt").T
x_p0, p0 = ler_dados("pressao_estag.txt").T
x_T, T = ler_dados("temperatura.txt").T
x_T0, T0 = ler_dados("temp_estag.txt").T
x_M, M = ler_dados("mach.txt").T
_, res1 = ler_dados("residuo1.txt").T
_, res2 = ler_dados("residuo2.txt").T
_, res3 = ler_dados("residuo3.txt").T

# Estetica do Grafico
fontSize = 14
fontSizeLegend = 12
darkBlue = (0.0, 0.129, 0.4784)
darkRed = (0.7176, 0.0705, 0.207)
plt.rc('font', family = 'serif')

# Eixo temporal
tempo = np.arange(len(res1))

# Plots
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Pressão
axs[0, 0].plot(x_p, p, color='blue')
axs[0, 0].set_title("Distribuicao de Pressão", fontsize=fontSize)
axs[0, 0].set_xlabel("$x$", fontsize=fontSize)
axs[0, 0].set_ylabel("$p$ (adimensional)", fontsize=fontSize)
axs[0, 0].grid(True)

# Pressão de Estagnação
axs[0, 1].plot(x_p0, p0, color='green')
axs[0, 1].set_title("Distribuicao de Pressão de Estagnação", fontsize=fontSize)
axs[0, 1].set_xlabel("$x$", fontsize=fontSize)
axs[0, 1].set_ylabel("$p_0$ (adimensional)", fontsize=fontSize)
axs[0, 1].set_ylim(0, 1.2)
axs[0, 1].grid(True)

# Temperatura
axs[1, 0].plot(x_T, T, color='red')
axs[1, 0].set_title("Distribuicao de Temperatura", fontsize=fontSize)
axs[1, 0].set_xlabel("$x$", fontsize=fontSize)
axs[1, 0].set_ylabel("$T$ (adimensional)", fontsize=fontSize)
axs[1, 0].grid(True)

# Mach
axs[1, 1].plot(x_M, M, color='purple')
axs[1, 1].set_title("Distribuicao de Número de Mach", fontsize=fontSize)
axs[1, 1].set_xlabel("$x$", fontsize=fontSize)
axs[1, 1].set_ylabel("$M$", fontsize=fontSize)
axs[1, 1].grid(True)

plt.tight_layout()
plt.show()

# Plots Separados

# Pressão
plt.figure(figsize=(8, 5))
plt.plot(x_p, p, color='darkBlue')
plt.title("Distribuicao de Pressão", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$p$ (adimensional)", fontsize=fontSize)
plt.grid(True)
plt.tight_layout()
plt.show()

# Pressão de Estagnação
plt.figure(figsize=(8, 5))
plt.plot(x_p0, p0, color='darkBlue')
plt.title("Distribuicao de Pressão de Estagnação", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$p_0$ (adimensional)", fontsize=fontSize)
plt.ylim(0, 1.2)
plt.grid(True)
plt.tight_layout()
plt.show()

# Temperatura
plt.figure(figsize=(8, 5))
plt.plot(x_T, T, color='darkBlue')
plt.title("Distribuicao de Temperatura", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$T$ (adimensional)", fontsize=fontSize)
plt.grid(True)
plt.tight_layout()
plt.show()

# Temperatura de Estagnação
plt.figure(figsize=(8, 5))
plt.plot(x_T0, T0, color='darkBlue')
plt.title("Distribuicao de Temperatura de Estagnação", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$T_0$ (adimensional)", fontsize=fontSize)
plt.ylim(0, 1.2)
plt.grid(True)
plt.tight_layout()
plt.show()

# Mach
plt.figure(figsize=(8, 5))
plt.plot(x_M, M, color='darkBlue')
plt.title("Distribuicao de Número de Mach", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$M$", fontsize=fontSize)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot dos Residuos
log_res1 = np.log10(res1)
log_res2 = np.log10(res2)
log_res3 = np.log10(res3)

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(tempo, log_res1, label='Resíduo - Massa Especifica')
ax.plot(tempo, log_res2, label='Resíduo - Quantidade de Movimento')
ax.plot(tempo, log_res3, label='Resíduo - Energia Total')

ax.set_title("Histórico de Convergência dos Resíduos", fontsize=fontSize)
ax.set_xlabel("Time Step", fontsize=fontSize)
ax.set_ylabel("$\log_{10}$(Resíduo)", fontsize=fontSize)
ax.legend(fontsize=fontSizeLegend)
ax.grid(True)

plt.tight_layout()
plt.show()
