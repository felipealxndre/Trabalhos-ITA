import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# NOVOS DADOS EXTRAÍDOS da Surface #1 (Asa Principal - Semi-Asa Positiva)
# Extraindo Yle (Posição na Semi-Asa [m])
Yle = np.array([
    0.0226, 0.2030, 0.5592, 1.0825, 1.7600, 2.5750, 3.5025, 4.5186, 5.6031, 6.7295,
    7.8699, 8.9962, 10.0808, 11.1925, 12.2923, 13.2536, 14.0526, 14.6133, 14.9501, 15.1206
])

# Extraindo cl_norm (Coeficiente de Sustentação Normal da Seção)
cl_norm = np.array([
    1.1776, 1.1738, 1.1926, 1.2358, 1.3035, 1.3741, 1.3980, 1.4142, 1.5142, 1.5643,
    1.5985, 1.6230, 1.6357, 1.6386, 1.6180, 1.5497, 1.3973, 1.1218, 0.7403, 0.2630
])


plt.figure(figsize=(10, 6))

# Encontrar o ponto de cl_norm máximo para marcar (simulando o ponto de estol)
cl_max_value = cl_norm.max()
idx_cl_max = cl_norm.argmax()
Yle_max = Yle[idx_cl_max]

plt.plot(Yle, cl_norm, marker='o', linestyle='-', color='black', label='Distribuição de $c_{l, norm}$ na Asa', markersize=3)


cl_max_ref = 1.64
plt.axhline(y=cl_max_ref, color='r', linestyle='--', label=f'$c_{{l, max}}$ (Referência da Figura) = {cl_max_ref}')

# Destacar o ponto de cl_norm máximo calculado nos dados
plt.plot(Yle_max, cl_max_value, marker='s', markersize=8, color='red', 
         label=f'Ponto de $c_{{l, norm}}$ máximo: {cl_max_value:.4f} @ Yle={Yle_max:.4f}m')

# Configurações do gráfico
plt.title('Distribuição do Coeficiente de Sustentação Normal $c_{l, norm}$ versus Posição na Semi-Asa')
plt.xlabel('Posição na Semi-Asa $Y_{le}$ [m]')
plt.ylabel('Coeficiente de Sustentação Normal da Seção $c_{l, norm}$')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.ylim(0.2, 1.8) 
plt.xlim(0, 16)
plt.show()