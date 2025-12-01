import design_tools
import numpy as np

# --- 1. CONFIGURAÇÃO INICIAL (Mesmos parâmetros do script anterior) ---
MTOW_GUESS = 434376.0
T0_GUESS = 100000.0
ALT_CRUISE = 35000 * 0.3048
MACH_CRUISE = 0.78
RANGE_CRUISE = 2000 * 1000
LOITER_TIME = 45 * 60
ALT_ALTCRUISE = 15000 * 0.3048
MACH_ALTCRUISE = 0.50
RANGE_ALTCRUISE = 200 * 1000

# --- 2. CARREGAR AERONAVE E GEOMETRIA ---
my_plane = design_tools.my_aircraft()
# Importante: Calcular geometria derivada (cordas médias, posições, etc)
new_dimensions = design_tools.geometry(my_plane)
my_plane['dimensions'].update(new_dimensions)

# --- 3. RODAR CÁLCULO DE PESO ---
W0, We, W_fuel, Mf_cruise, xcg_e, weightsvec = design_tools.weight(
    my_plane, MTOW_GUESS, T0_GUESS, ALT_CRUISE, MACH_CRUISE, RANGE_CRUISE,
    LOITER_TIME, ALT_ALTCRUISE, MACH_ALTCRUISE, RANGE_ALTCRUISE
)

# --- 4. EXTRAIR DADOS SOLICITADOS ---

# A) PESO DO MOTOR
# weightsvec = [W_w, W_h, W_v, W_f, W_nlg, W_mlg, W_eng_installed, W_allelse]
# O índice 6 é W_eng_installed (Peso total de todos motores instalados + naceles)
W_eng_total_installed = weightsvec[6]
n_engines = my_plane['data']['engines']['n']

# Para o diagrama de cortante, precisamos do peso de UM conjunto motor+nacele
W_eng_single = W_eng_total_installed / n_engines

# B) CENTROS AERODINÂMICOS (X_AC)
# O design_tools calcula Xac como: x_media + 0.25 * corda_media
# Vamos puxar esses dados calculados pela função geometry()
dim_w = my_plane['dimensions']['wing']
dim_h = my_plane['dimensions']['EH']

# X_AC da Asa
x_ac_w = dim_w['xm'] + 0.25 * dim_w['cm']

# X_AC da Empenagem Horizontal
x_ac_h = dim_h['xm'] + 0.25 * dim_h['cm']

# --- 5. EXIBIR RESULTADOS ---
print("-" * 40)
print("DADOS FALTANTES PARA ANÁLISE DE CARGAS")
print("-" * 40)
print(f"1. PESO DO MOTOR (W_eng)")
print(f"   - Peso Total Instalado ({int(n_engines)} motores): {W_eng_total_installed:.2f} N")
print(f"   - Peso por Motor (Para diagrama de Cortante): {W_eng_single:.2f} N")
print(f"     (Use este valor como a carga pontual em y = {my_plane['dimensions']['nacelle']['yn']} m)")
print("-" * 40)
print(f"2. CENTROS AERODINÂMICOS (X_ac)")
print(f"   - Asa (x_ac,w): {x_ac_w:.4f} m")
print(f"   - Empenagem Horizontal (x_ac,h): {x_ac_h:.4f} m")
print("-" * 40)
print("Copie estes valores para o seu documento ou memorial de cálculo.")