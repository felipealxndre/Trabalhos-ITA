import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d

# =============================================================================
# 1. DADOS DE ENTRADA E CONSTANTES
# =============================================================================
print("=== ANÁLISE DE CARGAS CORRIGIDA (V2) - PINGUINS DE MADAGASCAR ===\n")

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
sort_idx = np.argsort(y_vals)
y_vals = y_vals[sort_idx]
c_vals = c_vals[sort_idx]
c_cl_vals = c_cl_vals[sort_idx]
cm_c4_vals = cm_c4_vals[sort_idx]
area_vals = area_vals[sort_idx]

# =============================================================================
# 3. TRIMAGEM
# =============================================================================
M_ac_w_strip = cm_c4_vals * Q_DYN * area_vals * c_vals
M_ac_w_total = np.sum(M_ac_w_strip) * 2

lw = X_AC_W - X_CG
lh = X_AC_H - X_CG
A_mat = np.array([[1, 1], [lw, lh]])
B_vec = np.array([N_LOAD * W_TOTAL, -M_ac_w_total])
solution = np.linalg.solve(A_mat, B_vec)
L_W_trim = solution[0]

print(f"Sustentação Trimada Asa (L_W): {L_W_trim:.2f} N")

# =============================================================================
# 4. CÁLCULO DAS DISTRIBUIÇÕES BASE
# =============================================================================
# L'(y) Base (AVL) ajustada
L_prime_avl = Q_DYN * c_cl_vals
L_semi_avl_integrated = np.trapz(L_prime_avl, y_vals)
scale_factor = L_W_trim / (L_semi_avl_integrated * 2)
L_prime = L_prime_avl * scale_factor

# W'(y) Inercial
k_total = K_STR + K_FUEL
W_prime = N_LOAD * k_total * c_vals

# N'(y) Líquida (Distribuída)
N_prime = L_prime - W_prime

# =============================================================================
# 5. REFINAMENTO DA MALHA PARA PLOTAGEM (FIXAR O DEGRAU)
# =============================================================================
# Para criar o degrau vertical perfeito, inserimos o ponto Y_ENG duas vezes
# Criamos um vetor y aumentado para plotagem
y_augmented = np.sort(np.concatenate([y_vals, [Y_ENG, Y_ENG]]))

# Interpolamos a carga distribuída N' para esses novos pontos
# N' é contínua, então interpolação linear funciona bem
f_N = interp1d(y_vals, N_prime, kind='linear', fill_value="extrapolate")
N_aug = f_N(y_augmented)

# Integração da Cortante Distribuída (Sem Motor ainda)
# Integral de ponta a raiz (invertemos para integrar)
y_rev_aug = y_augmented[::-1]
N_rev_aug = N_aug[::-1]

# V_dist (apenas cargas distribuídas)
# Integral_{y}^{tip} N dy = - Integral_{tip}^{y} N dy
V_dist_raw = cumulative_trapezoid(N_rev_aug, y_rev_aug, initial=0)
V_dist_aug = -V_dist_raw[::-1]  # Inverte de volta para corresponder a y_augmented

# Adição da Carga Pontual do Motor com Lógica de Degrau
# O vetor y_augmented tem y crescente: [0, ..., 4.2, 4.2, ..., 16]
# O primeiro 4.2 é o limite esquerdo (inboard), o segundo é o direito (outboard).
# Fisicamente:
# Região Inboard (y < 4.2): V_total = V_dist - P_motor (tem a carga acumulada)
# Região Outboard (y > 4.2): V_total = V_dist (ainda não passou pelo motor vindo da ponta)
# No ponto exato 4.2:
#   Ponto "Esquerdo" (conectado à raiz) -> Sente o motor.
#   Ponto "Direito" (conectado à ponta) -> Não sente o motor.

P_eng_inertial = N_LOAD * W_ENG_UNIT
V_final_aug = np.zeros_like(V_dist_aug)

# Precisamos identificar qual dos dois pontos repetidos é qual.
# Como y_augmented está ordenado [..., 4.2, 4.2, ...], o índice menor é o da esquerda.
found_motor_start = False

for i, y in enumerate(y_augmented):
    val = V_dist_aug[i]

    # Lógica precisa para o degrau:
    if y < Y_ENG:
        # Claramente inboard
        val -= P_eng_inertial
    elif y == Y_ENG:
        # Estamos no ponto do motor.
        if not found_motor_start:
            # Primeira vez que vemos 4.2 (lado ESQUERDO/INBOARD do corte)
            # Deve incluir a carga do motor para conectar com a raiz
            val -= P_eng_inertial
            found_motor_start = True
        else:
            # Segunda vez que vemos 4.2 (lado DIREITO/OUTBOARD do corte)
            # Não inclui carga, conecta com a ponta
            pass
            # Se y > Y_ENG, não faz nada (mantém apenas carga distribuída)

    V_final_aug[i] = val

# Momento Fletor (Integrar a nova cortante refinada)
# M = Integral_{y}^{tip} V dy
V_final_rev = V_final_aug[::-1]
M_dist_raw = cumulative_trapezoid(V_final_rev, y_rev_aug, initial=0)
M_final_aug = -M_dist_raw[::-1]

# Interpolação para exibir as distribuições originais L' e W' no gráfico (apenas visual)
f_L = interp1d(y_vals, L_prime, kind='linear', fill_value="extrapolate")
L_aug = f_L(y_augmented)
f_W = interp1d(y_vals, W_prime, kind='linear', fill_value="extrapolate")
W_aug = f_W(y_augmented)

# =============================================================================
# 6. PLOTAGEM
# =============================================================================
print("\n--- GERANDO 3 GRÁFICOS SEPARADOS (V2 - DEGRAU CORRIGIDO) ---")

# Plot 1: Distribuições
plt.figure(1, figsize=(10, 6))
plt.plot(y_augmented, L_aug, label="Sustentação L'(y)", color='blue')
plt.plot(y_augmented, W_aug, label="Peso Inercial W'(y)", color='red', linestyle='--')
plt.plot(y_augmented, N_aug, label="Carga Líquida N'(y)", color='green', linewidth=2)
plt.title(f"Distribuição de Cargas (n={N_LOAD})")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Força Distribuída [N/m]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Plot 2: Força Cortante V(y)
plt.figure(2, figsize=(10, 6))
plt.plot(y_augmented, V_final_aug / 1000, color='purple', linewidth=2, label='Cortante Total')
plt.axvline(x=Y_ENG, color='orange', linestyle='--', alpha=0.5, label=f'Posição Motor ({Y_ENG}m)')
plt.title("Diagrama de Força Cortante V(y)")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Cortante [kN]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

# Plot 3: Momento Fletor M(y)
plt.figure(3, figsize=(10, 6))
plt.plot(y_augmented, M_final_aug / 1000, color='brown', linewidth=2, label='Momento Fletor')
plt.fill_between(y_augmented, M_final_aug / 1000, color='brown', alpha=0.1)
plt.title("Diagrama de Momento Fletor M(y)")
plt.xlabel("Semi-envergadura [m]")
plt.ylabel("Momento [kN.m]")
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()