import numpy as np
import matplotlib.pyplot as plt
import design_tools

# =============================================================================
# 1. DEFINIÇÃO DE CONSTANTES E PARÂMETROS
# =============================================================================
# Dados Físicos e de Voo
RHO = 1.225  # Densidade do ar ao nível do mar (kg/m³)
V_A_KT = 220.89  # Velocidade de manobra em nós
KT2MS = 0.514444  # Conversão de nós para m/s
V_A = V_A_KT * KT2MS  # Velocidade em m/s (~113.64 m/s)
Q = 0.5 * RHO * (V_A ** 2)  # Pressão dinâmica (Pa)
G = 9.81  # Aceleração da gravidade

# Parâmetros de Missão Estimados para o cálculo de peso
MTOW_GUESS = 434376.0
T0_GUESS = 100000.0
ALT_CRUISE_FT = 35000
MACH_CRUISE = 0.78
RANGE_CRUISE_M = 2000 * 1000
LOITER_TIME_S = 45 * 60
ALT_ALTCRUISE_FT = 15000
MACH_ALTCRUISE = 0.50
RANGE_ALTCRUISE_M = 200 * 1000

# =============================================================================
# 2. CÁLCULO DAS MASSAS (Usando design_tools.py)
# =============================================================================
print("--- CALCULANDO MASSAS USANDO DESIGN_TOOLS ---")

# Carrega a aeronave
my_plane = design_tools.my_aircraft()

# --- CORREÇÃO: CALCULAR GEOMETRIA ANTES DE CHAMAR WEIGHT ---
new_dimensions = design_tools.geometry(my_plane)
my_plane['dimensions'].update(new_dimensions)
# -----------------------------------------------------------

# Executa a iteração de peso
W0, We, W_fuel, Mf_cruise, xcg_e, weightsvec = design_tools.weight(
    my_plane,
    MTOW_GUESS,
    T0_GUESS,
    ALT_CRUISE_FT * 0.3048,
    MACH_CRUISE,
    RANGE_CRUISE_M,
    LOITER_TIME_S,
    ALT_ALTCRUISE_FT * 0.3048,
    MACH_ALTCRUISE,
    RANGE_ALTCRUISE_M
)

W_wing = weightsvec[0]
S_wing = my_plane['geo_param']['wing']['S']

# Cálculo das constantes de distribuição trapezoidal
k_structure = W_wing / S_wing
k_fuel = W_fuel / S_wing

print(f"MTOW Calculado (W0): {W0:.2f} N")
print(f"Peso da Asa (Estrutural): {W_wing:.2f} N")
print(f"Peso do Combustível: {W_fuel:.2f} N")
print("-" * 30)
print(f"Fator k_estrutura (W_wing/S): {k_structure:.2f} N/m²")
print(f"Fator k_combustível (W_fuel/S): {k_fuel:.2f} N/m²")
print("=" * 60)

# =============================================================================
# 3. LEITURA DO ARQUIVO FORCES.TXT E PLOTAGEM
# =============================================================================
print("\n--- PROCESSANDO ARQUIVO FORCES.TXT ---")

file_path = 'forces.txt'

y_vals = []
c_cl_vals = []
cl_vals = []  # Lista para armazenar o Cl local

reading_surface_1 = False

try:
    with open(file_path, 'r') as f:
        for line in f:
            clean_line = line.strip()

            # Identifica o início da Superfície 1 (Asa Principal)
            if "Surface # 1" in clean_line and "Asa Principal" in clean_line:
                reading_surface_1 = True
                continue

            # Se encontrar outra superfície, para de ler
            if "Surface #" in clean_line and reading_surface_1:
                if "Surface # 1" not in clean_line:
                    break

            # Ler os dados da tabela
            if reading_surface_1:
                try:
                    parts = clean_line.split()
                    if len(parts) > 10 and parts[0].isdigit():
                        # Colunas AVL forces.txt:
                        # 0: j, 1: Yle, 4: c cl, 7: cl
                        y = float(parts[1])
                        c_cl = float(parts[4])
                        cl = float(parts[7])  # <--- Leitura do Cl

                        y_vals.append(y)
                        c_cl_vals.append(c_cl)
                        cl_vals.append(cl)
                except ValueError:
                    continue

    if not y_vals:
        print("ERRO: Nenhum dado extraído do forces.txt.")
    else:
        # Converter para arrays numpy
        y_arr = np.array(y_vals)
        c_cl_arr = np.array(c_cl_vals)
        cl_arr = np.array(cl_vals)

        # Calcular Distribuição de Sustentação L'(y) [N/m]
        L_prime = Q * c_cl_arr

        # --- PLOT 1: Distribuição de Sustentação L'(y) ---
        plt.figure(figsize=(10, 6))
        plt.plot(y_arr, L_prime, label=f"Sustentação L'(y)", color='blue', linewidth=2)
        plt.fill_between(y_arr, L_prime, alpha=0.3, color='blue')

        plt.title('Distribuição de Sustentação na Asa (L\')', fontsize=14)
        plt.xlabel('Envergadura y [m]', fontsize=12)
        plt.ylabel("Sustentação L'(y) [N/m]", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        # Integral para verificação
        L_total_calc = np.trapz(L_prime, y_arr)
        plt.text(0.05, 0.95, f'Sustentação Integrada (Semi-asa): {L_total_calc / 1000:.1f} kN',
                 transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))

        # --- PLOT 2: Distribuição de Cl ---
        plt.figure(figsize=(10, 6))
        plt.plot(y_arr, cl_arr, label="Cl local", color='red', linewidth=2, linestyle='-')

        plt.title('Distribuição de Coeficiente de Sustentação (Cl) na Asa', fontsize=14)
        plt.xlabel('Envergadura y [m]', fontsize=12)
        plt.ylabel("Coeficiente de Sustentação Cl [-]", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        plt.show()
        print("Gráficos gerados com sucesso.")

except FileNotFoundError:
    print(f"ERRO: Arquivo '{file_path}' não encontrado no diretório local.")