import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# =============================================================================
# 1. DADOS DE ENTRADA E CONSTANTES
# =============================================================================
print("=== ANÁLISE DE CARGAS CORRIGIDA - PINGUINS DE MADAGASCAR ===\n")

# Parâmetros Físicos e de Voo
N_LOAD = 2.5
RHO = 1.225
V_A = 113.64
Q_DYN = 0.5 * RHO * V_A ** 2

# Pesos (Newtons)
W_TOTAL = 410443.36
W_WING = 41036.98
W_FUEL = 71370.75
W_ENG_UNIT = 11068.29
N_ENGINES = 2

# Geometria
S_REF = 105.5
X_CG = 14.567
X_AC_W = 15.6933
X_AC_H = 28.3901
Y_ENG = 4.2

# Fatores de Distribuição
K_STR = W_WING / S_REF
K_FUEL = W_FUEL / S_REF


# =============================================================================
# 2. PROCESSAMENTO DO ARQUIVO FORCES.TXT
# =============================================================================
def read_avl_forces(filename):
    y_list, chord_list, c_cl_list, cm_c4_list, area_list = [], [], [], [], []
    reading = False
    try:
        with open(filename, 'r') as f:
            for line in f:
                if "Surface # 1" in line and "Asa Principal" in line:
                    reading = True
                    continue
                if "Surface #" in line and reading and "Surface # 1" not in line:
                    break

                if reading:
                    try:
                        parts = line.split()
                        if len(parts) > 10 and parts[0].isdigit():
                            y_list.append(float(parts[1]))
                            chord_list.append(float(parts[2]))
                            area_list.append(float(parts[3]))
                            c_cl_list.append(float(parts[4]))
                            cm_c4_list.append(float(parts[10]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"ERRO: Arquivo '{filename}' não encontrado.")
        exit()

    return np.array(y_list), np.array(chord_list), np.array(c_cl_list), np.array(cm_c4_list), np.array(area_list)


y_vals, c_vals, c_cl_vals, cm_c4_vals, area_vals = read_avl_forces('forces.txt')

# Ordenar arrays da Raiz (y=0) para a Ponta (y=b/2)
# Isso garante que o gráfico comece em 0 na esquerda
sort_idx = np.argsort(y_vals)
y_vals = y_vals[sort_idx]
c_vals = c_vals[sort_idx]
c_cl_vals = c_cl_vals[sort_idx]
cm_c4_vals = cm_c4_vals[sort_idx]
area_vals = area_vals[sort_idx]

# =============================================================================
# 3. TRIMAGEM
# =============================================================================
# Momento Aerodinâmico (Asa inteira)
M_ac_w_strip = cm_c4_vals * Q_DYN * area_vals * c_vals
M_ac_w_total = np.sum(M_ac_w_strip) * 2

# Sistema Linear
lw = X_AC_W - X_CG
lh = X_AC_H - X_CG
A_mat = np.array([[1, 1], [lw, lh]])
B_vec = np.array([N_LOAD * W_TOTAL, -M_ac_w_total])
solution = np.linalg.solve(A_mat, B_vec)
L_W_trim = solution[0]

print(f"Sustentação Trimada Asa (L_W): {L_W_trim:.2f} N")

# =============================================================================
# 4. DISTRIBUIÇÃO DE CARGAS
# =============================================================================
# L'(y) Base (AVL)
L_prime_avl = Q_DYN * c_cl_vals

# Escalar para bater com Trimagem
L_semi_avl_integrated = np.trapz(L_prime_avl, y_vals)
L_total_avl_integrated = L_semi_avl_integrated * 2
scale_factor = L_W_trim / L_total_avl_integrated
L_prime = L_prime_avl * scale_factor

# W'(y) Inercial
k_total = K_STR + K_FUEL
W_prime = N_LOAD * k_total * c_vals

# N'(y) Líquida
N_prime = L_prime - W_prime

# =============================================================================
# 5. INTEGRAÇÃO (DIAGRAMAS V e M)
# =============================================================================
# Para integrar da Ponta para Raiz, invertemos os arrays
y_rev = y_vals[::-1]  # Vai de ~16.4 até 0
N_prime_rev = N_prime[::-1]

# Cortante V(y)
# Integral de N' dy. Como y_rev decresce, o dy implícito é negativo.
# cumulative_trapezoid(N, y_rev) calcula Integral_{16}^{y} N dy = - Integral_{y}^{16} N dy
# O roteiro pede V(y) = Integral_{y}^{16} N dy.
# Logo, V_dist = - cumulative_trapezoid(...)
V_dist_raw = cumulative_trapezoid(N_prime_rev, y_rev, initial=0)
V_dist = -V_dist_raw

# Adicionar Motor (Carga Pontual para BAIXO em y=4.2)
P_eng_inertial = N_LOAD * W_ENG_UNIT
V_total = np.zeros_like(V_dist)

for i, y in enumerate(y_rev):
    v_val = V_dist[i]
    # Se estamos à esquerda do motor (mais perto da raiz), somamos a carga do motor.
    # Motor aponta para baixo (-). Cortante acumula cargas.
    if y < Y_ENG:
        v_val -= P_eng_inertial
    V_total[i] = v_val

# Momento M(y)
# M(y) = Integral_{y}^{16} V(eta) deta
# Novamente, integrando sobre y_rev (decrescente), resultado sai negativo da integral desejada.
M_dist_raw = cumulative_trapezoid(V_total, y_rev, initial=0)
M_dist = -M_dist_raw

# =============================================================================
# 6. PLOTAGEM CORRIGIDA
# =============================================================================
print("\n--- GERANDO 3 GRÁFICOS SEPARADOS ---")

# PREPARAÇÃO DOS DADOS PARA PLOTAGEM (EIXO X: 0 -> 16)
# y_vals já é 0 -> 16.
# L_prime, W_prime, N_prime já estão alinhados com y_vals (0->16).

# V_total e M_dist foram calculados em y_rev (16 -> 0).
# Precisamos inverter eles de volta para plotar contra y_vals (0 -> 16).
V_plot = V_total[::-1]
M_plot = M_dist[::-1]

# Plot 1: Distribuições
plt.figure(1, figsize=(10, 6))
plt.plot(y_vals, L_prime, label="Sustentação L'(y)", color='blue')
plt.plot(y_vals, W_prime, label="Peso Inercial W'(y)", color='red', linestyle='--')
plt.plot(y_vals, N_prime, label="Carga Líquida N'(y)", color='green', linewidth=2)
plt.title(f"Distribuição de Cargas (n={N_LOAD})")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Força Distribuída [N/m]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Plot 2: Força Cortante V(y)
plt.figure(2, figsize=(10, 6))
plt.plot(y_vals, V_plot / 1000, color='purple', linewidth=2, label='Cortante Total')
# Linha vertical do motor apenas para referência visual (sem seta)
plt.axvline(x=Y_ENG, color='orange', linestyle='--', alpha=0.5, label=f'Posição Motor ({Y_ENG}m)')
plt.title("Diagrama de Força Cortante V(y)")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Cortante [kN]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Plot 3: Momento Fletor M(y)
plt.figure(3, figsize=(10, 6))
plt.plot(y_vals, M_plot / 1000, color='brown', linewidth=2, label='Momento Fletor')
plt.fill_between(y_vals, M_plot / 1000, color='brown', alpha=0.1)
plt.title("Diagrama de Momento Fletor M(y)")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Momento [kN.m]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()