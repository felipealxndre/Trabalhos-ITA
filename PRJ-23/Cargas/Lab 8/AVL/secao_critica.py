import pandas as pd
import matplotlib.pyplot as plt
import io
import numpy as np

# Dados fornecidos
data = """
Yle cl_norm cl
0.0245 0.6893 0.5836
0.2199 0.6851 0.5800
0.6058 0.6991 0.5919
1.1727 0.7298 0.6179
1.9067 0.7767 0.6576
2.7897 0.8191 0.6935
3.7943 0.8181 0.6927
4.8950 0.8588 0.7271
6.0698 0.9115 0.7717
7.2899 0.9433 0.7986
8.5252 0.9643 0.8164
9.7452 0.9823 0.8317
10.9201 0.9934 0.8411
12.1242 1.0000 0.8467
13.3154 0.9962 0.8435
14.3565 0.9654 0.8174
15.2219 0.8830 0.7476
15.8292 0.7183 0.6081
16.1941 0.4754 0.4025
16.3788 0.1683 0.1425
"""

# Ler os dados
df = pd.read_csv(io.StringIO(data), sep='\s+')

# Ponto de estol (14º ponto, índice 13 em Python)
stall_index = 13
stall_yle = df['Yle'].iloc[stall_index]
stall_cl = df['cl'].iloc[stall_index]
stall_alpha = 8.676

# Valor constante de clmax do aerofólio
cl_max_value = 0.8467

# Configuração do gráfico
plt.figure(figsize=(10, 6))

# 1. Curva principal: cl vs Yle em preto
plt.plot(df['Yle'], df['cl'], color='black', marker='o', linestyle='-', linewidth=2, label=r'Distribuição de $c_l$ na Asa')

# 2. Reta horizontal de clmax em vermelho
plt.plot(df['Yle'], np.full_like(df['Yle'], cl_max_value), color='red', linestyle='--', linewidth=1.5, label=r'$c_{l_{max}} = ' + str(cl_max_value) + r'$')

# 3. Ponto focal de estol (índice 13)
plt.plot(stall_yle, stall_cl, marker='s', markersize=10, color='red', markeredgecolor='black', linewidth=0,
         label=r'Ponto de estol, $\alpha = ' + str(stall_alpha) + r'^\circ$')


# Configurações de rótulos e título
plt.title(r'Distribuição do Coeficiente de Sustentação $c_l$ versus Posição na Semi-Asa', fontsize=14)
plt.xlabel(r'Posição na Semi-Asa $Y_{le}$ [m]', fontsize=12)
plt.ylabel(r'Coeficiente de Sustentação do Aerofólio $c_l$', fontsize=12)

# Adicionar grade e legenda
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10, loc='lower left')
plt.show()