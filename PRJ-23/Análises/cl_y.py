import matplotlib.pyplot as plt
import numpy as np

# Dados da Semi-Asa Positiva (Surface #1: Strips 1 a 20)
# Yle (Posição na Semi-Asa [m])
Yle = np.array([
    0.0200, 0.1788, 0.4927, 0.9538, 1.5507, 2.2689, 3.0861, 3.9814, 4.9370, 5.9295,
    6.9343, 7.9267, 8.8824, 9.8618, 10.8306, 11.6773, 12.3812, 12.8752, 13.1719, 13.3221
])

# cl (Coeficiente de Sustentação da Seção)
cl = np.array([
    1.1187, 1.1311, 1.1473, 1.1861, 1.2504, 1.3268, 1.3776, 1.3603, 1.4558, 1.5354,
    1.5829, 1.6150, 1.6331, 1.6400, 1.6223, 1.5532, 1.3984, 1.1204, 0.7384, 0.2651
])


plt.figure(figsize=(10, 6))

# Encontrar o ponto de cl máximo para marcar (simulando o ponto de estol)
cl_max_value = cl.max()
idx_cl_max = cl.argmax()
Yle_max = Yle[idx_cl_max]

# Plotar a distribuição do cl
plt.plot(Yle, cl, marker='o', linestyle='-', color='black', label='Distribuição de $c_l$ na Asa', markersize=3)

# Adicionar a linha horizontal do cl_max de referência da figura (1.64)
cl_max_ref = 1.64
plt.axhline(y=cl_max_ref, color='r', linestyle='--', label=f'$c_{{l, max}}$ = {cl_max_ref}')

# Destacar o ponto de cl máximo calculado nos dados
plt.plot(Yle_max, cl_max_value, marker='s', markersize=8, color='red', label=f'Ponto de $c_l$ máximo: {cl_max_value:.4f}')

# Configurações do gráfico
plt.title('Distribuição do Coeficiente de Sustentação $c_l$ versus Posição na Semi-Asa')
plt.xlabel('Posição na Semi-Asa $Y_{le}$ [m]')
plt.ylabel('Coeficiente de Sustentação da Seção $c_l$')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.ylim(0.2, 1.8) # Limites para replicar o visual da sua figura
plt.show()
