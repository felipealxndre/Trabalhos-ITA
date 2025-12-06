import math
import numpy as np
import matplotlib.pyplot as plt
import design_tools as dt
from aux_tools import atmosphere
from scipy.integrate import cumulative_trapezoid
import os

g = 9.81

# =============================================================================
# Configuração da aeronave
# =============================================================================
aircraft = dt.my_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)

S = aircraft['geo_param']['wing']['S']
c_mac = aircraft['dimensions']['wing']['cm']
b_span = aircraft['dimensions']['wing']['b']
taper_ratio = aircraft['geo_param']['wing']['taper']

x_ac_w = aircraft['dimensions']['wing']['xm'] + 0.25 * aircraft['dimensions']['wing']['cm']
x_ac_h = aircraft['dimensions']['EH']['xm'] + 0.25 * aircraft['dimensions']['EH']['cm']

rho_SL = 1.2250
rho_11000 = 0.9217

def compute_Cm_acw_from_avl(filename, S_ref, c_mac):
    y_list, chord_list, area_list, cm_c4_list = [], [], [], []
    reading = False

    possible_paths = [
        filename,
        os.path.join('Cargas', filename),
        os.path.join(os.path.dirname(__file__), filename),
        os.path.join(os.path.dirname(__file__), '..', filename),
    ]

    file_found = None
    for path in possible_paths:
        if os.path.exists(path):
            file_found = path
            break

    if file_found is None:
        raise FileNotFoundError(f"{filename} não encontrado para cálculo de Cm_acw.")

    with open(file_found, 'r') as f:
        for line in f:
            if "Surface # 1" in line and "Asa Principal" in line:
                reading = True
                continue
            if "Surface #" in line and reading and "Surface # 1" not in line:
                break

            if reading:
                parts = line.split()
                if len(parts) > 10 and parts[0].isdigit():
                    y_list.append(float(parts[1]))
                    chord_list.append(float(parts[2]))
                    area_list.append(float(parts[3]))
                    cm_c4_list.append(float(parts[10]))

    chord_arr = np.array(chord_list)
    area_arr = np.array(area_list)
    cm_c4_arr = np.array(cm_c4_list)

    # Momento total no AC da asa (asa inteira)
    M_ac_w_non_dim = 2.0 * np.sum(cm_c4_arr * area_arr * chord_arr)

    # Coeficiente global Cm_acw
    Cm_acw = M_ac_w_non_dim / (S_ref * c_mac)
    return Cm_acw

Cm_acw = compute_Cm_acw_from_avl("forces_anv_trimada.txt", S, c_mac)
print(f"Cm_acw calculado a partir do AVL: {Cm_acw:.4f}")

W0_guess = 43090 * g
T0_guess = 125600

Mach_cruise = 0.75
altitude_cruise = 11000
range_cruise = 2_390_000.0

Mach_altcruise = 0.4
range_altcruise = 370_000.0
altitude_altcruise = 4572

loiter_time = 2700

pesos = dt.analyze(
    aircraft,
    W0_guess, T0_guess,
    Mach_cruise, altitude_cruise, range_cruise,
    Mach_altcruise, range_altcruise, altitude_altcruise,
    loiter_time,
    0, 1520, 0.34906585039887, 0,
    0, 1520, 0.69813170079773, 0,
    0.84
)

MTOW = pesos[0]
MZFW = (pesos[0] - pesos[1])

# Rodamos o weight isolado só para extrair o vetor de pesos estruturais
W0_weight, We_weight, Wf_weight, Mf_cruise_w, xcg_e_weight, weightsvec = dt.weight(
    aircraft,
    MTOW,              # usa o MTOW já convergido
    T0_guess,          # chuta o T0; o weight converge W0 internamente
    altitude_cruise, Mach_cruise, range_cruise,
    loiter_time,
    altitude_altcruise, Mach_altcruise, range_altcruise
)

x_cg = aircraft['weights']['xcg_fwd']

# Peso estrutural da asa (W_w)
W_wing_struct = weightsvec[0]

# Peso total de motores instalados (todos os motores + naceles)
W_eng_total_installed = weightsvec[6]

# Número de motores
n_engines = aircraft['data']['engines']['n']

# Peso de UM motor instalado (para usar como carga pontual por lado)
W_eng_single = W_eng_total_installed / n_engines

# Peso total de combustível (vem do analyze)
W_fuel_total = pesos[1]

# Fração asa+combustível em relação ao MTOW
FRAC_WING_FUEL = (W_wing_struct + W_fuel_total) / MTOW

# Fração de UM motor em relação ao MTOW
FRAC_ENG = W_eng_single / MTOW

print(f"FRAC_WING_FUEL = {FRAC_WING_FUEL:.5f}")
print(f"FRAC_ENG       = {FRAC_ENG:.5f}")

# Regional TBP's 
# Roskam
W_frac_engine_start = 0.990
W_frac_taxi = 0.990  
W_frac_take_off = 0.995 
W_frac_climb = 0.98 

frac_pre_cruise = (
    W_frac_engine_start *
    W_frac_taxi *
    W_frac_take_off *
    W_frac_climb 
)

W_cruise = MTOW * frac_pre_cruise

MLW = 0.84 * MTOW

# =============================================================================
# Geometria da corda local
# =============================================================================
def c_local(y):
    c_root = 2 * S / (b_span * (1 + taper_ratio))
    return c_root * (1 - (1 - taper_ratio) * (2 * y / b_span))

# =============================================================================
# Leitura AVL: c*cl em strips
# =============================================================================
def read_avl_forces(filename):
    y_list = []
    ccl_list = []
    reading = False

    possible_paths = [
        filename,
        os.path.join('Cargas', filename),
        os.path.join(os.path.dirname(__file__), filename),
        os.path.join(os.path.dirname(__file__), '..', filename),
    ]

    file_found = None
    for path in possible_paths:
        if os.path.exists(path):
            file_found = path
            break

    if file_found is None:
        raise FileNotFoundError("forces_anv_trimada.txt não encontrado.")

    with open(file_found, 'r') as f:
        for line in f:
            if "Surface # 1" in line and "Asa Principal" in line:
                reading = True
                continue
            if "Surface #" in line and reading and "Surface # 1" not in line:
                break

            if reading:
                parts = line.split()
                if len(parts) > 10 and parts[0].isdigit():
                    y_list.append(float(parts[1]))
                    ccl_list.append(float(parts[4]))

    y_arr = np.array(y_list)
    ccl_arr = np.array(ccl_list)

    sort_idx = np.argsort(y_arr)
    y_arr = y_arr[sort_idx]
    ccl_arr = ccl_arr[sort_idx]

    if y_arr[0] > 1e-3:
        y_arr = np.concatenate([[0.0], y_arr])
        ccl_arr = np.concatenate([[ccl_arr[0]], ccl_arr])

    return y_arr, ccl_arr


y_avl, ccl_avl = read_avl_forces("forces_anv_trimada.txt")

# =============================================================================
# Função de cálculo de cargas para UM caso
# =============================================================================
def calcular_caso(W_case, n_case, V_case, rho_case, nome):
    P_limit = n_case * W_case
    P_ult = 1.5 * P_limit

    q = 0.5 * rho_case * V_case**2
    M_ac_w = q * S * c_mac * Cm_acw

    lw = x_ac_w - x_cg
    lh = x_ac_h - x_cg

    A = np.array([[1, 1], [lw, lh]])
    b = np.array([P_limit, -M_ac_w])
    L_w, L_h = np.linalg.solve(A, b)

    y = y_avl
    ccl = ccl_avl

    L_semi = L_w / 2
    L_prime_base = q * ccl
    L_semi_base = np.trapz(L_prime_base, y)
    scale = L_semi / L_semi_base
    L_prime = L_prime_base * scale

    c_y = c_local(y)
    W_wing_fuel_total = FRAC_WING_FUEL * W_case
    W_wing_fuel_semi = W_wing_fuel_total / 2

    k_w = W_wing_fuel_semi / np.trapz(c_y, y)
    W_prime = k_w * c_y

    N_prime = L_prime - W_prime

    y_eng = aircraft['dimensions']['nacelle']['yn']
    W_eng = FRAC_ENG * W_case

    y_aug = np.sort(np.concatenate([y, [y_eng, y_eng]]))
    N_aug = np.interp(y_aug, y, N_prime)

    y_rev = y_aug[::-1]
    N_rev = N_aug[::-1]

    V_dist_raw = cumulative_trapezoid(N_rev, y_rev, initial=0)
    V_dist = -V_dist_raw[::-1]

    V_final = np.zeros_like(V_dist)
    found_start = False

    for i, yi in enumerate(y_aug):
        val = V_dist[i]
        if yi < y_eng:
            val += W_eng
        elif yi == y_eng and not found_start:
            val += W_eng
            found_start = True
        V_final[i] = val

    V_rev = V_final[::-1]
    M_raw = cumulative_trapezoid(V_rev, y_rev, initial=0)
    M_final = -M_raw[::-1]
    
    return {
        'y_aug': y_aug,
        'L_prime': L_prime,
        'W_prime': W_prime,
        'N_prime': N_prime,
        'N_aug': N_aug,
        'V_final': V_final,
        'M_final': M_final,
        'L_w': L_w,
        'L_h': L_h
    }


# =============================================================================
# Casos
# =============================================================================

casos = [
    ("MTOW @ Nível do Mar", MTOW, 2.5, 108.88, rho_SL),
    ("MZFW @ Nível do Mar", MZFW, 2.854, 164.62, rho_SL),
    ("Peso de Cruzeiro @ 11000 ft", W_cruise, 2.574, 164.62, rho_11000),
    ("MLW @ 11000 ft", MLW, 2.845, 164.62, rho_11000),
]

# Calcular e armazenar todos os resultados
resultados_casos = {}
for nome, Wc, nc, Vc, rhoc in casos:
    resultados = calcular_caso(Wc, nc, Vc, rhoc, nome)
    resultados_casos[nome] = {
        'W_case': Wc,
        'n_case': nc,
        'L_w': resultados['L_w'],
        'L_h': resultados['L_h'],
        'resultados': resultados
    }

# Imprimir resumo
print("\nCargas de balanceamento e incremento de carga na asa (valores em N):")
for nome, res in resultados_casos.items():
    Wc = res["W_case"]
    nc = res["n_case"]
    L_w = res["L_w"]
    L_h = res["L_h"]

    incremento_asa = L_w - nc * Wc
    n_asa = L_w / Wc

    print(f"{nome}: n = {nc:.3f}")
    print(f"  LH (empenagem horizontal) = {L_h:10.1f} N")
    print(f"  LW (asa)                  = {L_w:10.1f} N")
    print(f"  Incremento na asa LW - nW = {incremento_asa:10.1f} N")
    print(f"  Fator de carga na asa     = {n_asa:6.3f} g\n")

# Configuração de fontes para os gráficos
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11
})

# Gráficos
plt.figure(figsize=(12, 8))
for nome, res in resultados_casos.items():
    r = res['resultados']
    plt.plot(r['y_aug'], np.interp(r['y_aug'], y_avl, r['L_prime']), label=f"L'(y) - {nome}")
    plt.plot(r['y_aug'], np.interp(r['y_aug'], y_avl, r['W_prime']), label=f"W'(y) - {nome}", linestyle='-.')

plt.title("Distribuição de Cargas L(y), W(y) - Comparação")
plt.xlabel("y [m]")
plt.ylabel("Força [N/m]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('Cargas/distribuicao_cargas.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(12, 8))
for nome, res in resultados_casos.items():
    r = res['resultados']
    plt.plot(r['y_aug'], r['V_final']/1000, label=f"V(y) - {nome}")

plt.title("Força Cortante V(y) - Comparação")
plt.xlabel("y [m]")
plt.ylabel("V(y) [kN]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('Cargas/forca_cortante.png', dpi=300, bbox_inches='tight')

plt.figure(figsize=(12, 8))
for nome, res in resultados_casos.items():
    r = res['resultados']
    plt.plot(r['y_aug'], r['M_final']/1000, label=f"M(y) - {nome}")

plt.title("Momento Fletor M(y) - Comparação")
plt.xlabel("y [m]")
plt.ylabel("M(y) [kN·m]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('Cargas/momento_fletor.png', dpi=300, bbox_inches='tight')

plt.show()

