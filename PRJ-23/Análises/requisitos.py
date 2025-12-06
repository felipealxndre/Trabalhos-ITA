import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import design_tools_optmized as dt
import aux_tools as at
from plot3d import plot3d

# -----------------------------------------------------------------------------
# Definição da aeronave e geometria
# -----------------------------------------------------------------------------
aircraft = dt.my_aircraft()

new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)

gravity = 9.81

# -----------------------------------------------------------------------------
# Parâmetros de missão e certificação
# -----------------------------------------------------------------------------
T0_guess = 125600.0        # chute inicial de empuxo total (N)
W0_guess = 490000.0        # chute inicial de peso máximo de decolagem (N)

altitude_cruise   = 11000.0   # m
Mach_cruise       = 0.75
range_cruise      = 3700e3    # m (3700 km)

range_altcruise   = 370400.0  # m (200 NM)
altitude_altcruise = 4572.0   # m
Mach_altcruise    = 0.40

loiter_time       = 2700.0    # s (45 min)

# Requisitos de pista (entrada do problema)
distance_takeoff_req  = 1800.0   # m (requisito)
distance_landing_req  = 1150.0   # m (requisito)

# Deflexões de hipersustentadores
TO_flap_def = 20.0 * np.pi / 180.0
LD_flap_def = aircraft['data']['flap']['max_def']
TO_slat_def = 0.0
LD_slat_def = 0.0

# Outros parâmetros
h_ground        = 10.668  # m (~35 ft)
altitude_takeoff = 0.0    # ISA, SL
altitude_landing = 0.0    # ISA, SL
MLW_frac        = 0.84    # MLW/MTOW

# -----------------------------------------------------------------------------
# Análise principal (peso, desempenho, etc.)
# -----------------------------------------------------------------------------
W0, Wf, T0, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, \
frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = \
    dt.analyze(
        aircraft, W0_guess, T0_guess,
        Mach_cruise, altitude_cruise, range_cruise,
        Mach_altcruise, range_altcruise, altitude_altcruise,
        loiter_time,
        altitude_takeoff, distance_takeoff_req, TO_flap_def, TO_slat_def,
        altitude_landing, distance_landing_req, LD_flap_def, LD_slat_def,
        MLW_frac
    )

aircraft['weights']['W0'] = W0

T0_available, T0vec_req, S_wlan_req = dt.performance(
    aircraft,
    TO_flap_def, LD_flap_def,
    TO_slat_def, LD_slat_def,
    h_ground,
    altitude_takeoff, distance_takeoff_req,
    altitude_landing, distance_landing_req, MLW_frac,
    altitude_cruise, Mach_cruise
)

T0vec_req = np.array(T0vec_req)

# -----------------------------------------------------------------------------
# Cálculo explícito das distâncias de decolagem e pouso
# -----------------------------------------------------------------------------
S_w = aircraft['geo_param']['wing']['S']

# Atmosfera para decolagem (ISA, SL)
T_to, p_to, rho_to, mi_to = at.atmosphere(altitude_takeoff, 288.15)

# CLmax em decolagem (mesma condição usada em dt.performance)
_, CLmaxTO = dt.aerodynamics(
    aircraft,
    Mach=0.2,
    altitude=altitude_takeoff,
    n_engines_failed=0,
    flap_def=TO_flap_def,
    slat_def=TO_slat_def,
    lg_down=0,
    h_ground=h_ground,
    W0_guess=W0
)

sigma_to = rho_to / 1.225

# Fórmula usada em performance rearranjada para S_TO:
# T0/W0 = 0.2387 / (sigma * CLmaxTO * S_TO) * (W0 / S_w)
#  ⇒ S_TO = 0.2387 / (sigma * CLmaxTO) * (W0**2) / (S_w * T0_available)
distance_takeoff_calc = 0.2387 / (sigma_to * CLmaxTO) * (W0**2) / (S_w * T0_available)

# Atmosfera para pouso (ISA, SL)
T_ld, p_ld, rho_ld, mi_ld = at.atmosphere(altitude_landing, 288.15)

# CLmax em pouso (mesma condição usada em dt.performance)
_, CLmaxLD = dt.aerodynamics(
    aircraft,
    Mach=0.2,
    altitude=altitude_landing,
    n_engines_failed=0,
    flap_def=LD_flap_def,
    slat_def=LD_slat_def,
    lg_down=0,
    h_ground=h_ground,
    W0_guess=W0 * MLW_frac
)

# Coeficientes de Torenbeek (mesmos de dt.performance)
h_land = 15.3
f_land = 5.0 / 3.0
a_g    = 0.5

x_land  = 1.52 / a_g + 1.69
A_land  = gravity / f_land / x_land
B_land  = -10.0 * gravity * h_land / x_land

# Em dt.performance:
# S_wlan = W0*MLW_frac / (rho * (A_land * S_L + B_land) * CLmaxLD)
# Rearranjando para S_L (distância de pouso) com S_w conhecido:
# A_land * S_L + B_land = W0*MLW_frac / (rho * CLmaxLD * S_w)
# S_L = (W0*MLW_frac / (rho * CLmaxLD * S_w) - B_land) / A_land
distance_landing_calc = (
    (W0 * MLW_frac) / (rho_ld * CLmaxLD * S_w) - B_land
) / A_land

# -----------------------------------------------------------------------------
# Gradientes de subida e demais requisitos
# -----------------------------------------------------------------------------
b_wing = aircraft['dimensions']['wing']['b']
n_pax = 90
n_tripulantes = 5

excess_thrust = T0_available - T0vec_req
excess_grad   = excess_thrust / W0

grad_25_111 = 1.2 + (excess_grad[2] * 100.0)
grad_121_a  = 0.0 + (excess_grad[3] * 100.0)
grad_121_b  = 2.4 + (excess_grad[4] * 100.0)
grad_121_c  = 1.2 + (excess_grad[5] * 100.0)
grad_25_119 = 3.2 + ((excess_thrust[6] / (W0 * MLW_frac)) * 100.0)
grad_121_d  = 2.1 + ((excess_thrust[7] / (W0 * MLW_frac)) * 100.0)

excess_thrust_cruise = T0_available - T0vec_req[1]

# -----------------------------------------------------------------------------
# Matriz de conformidade
# -----------------------------------------------------------------------------
data = []

data.append({"Parâmetro": "Número de passageiros",
             "Requisito": "90",
             "Valor": n_pax,
             "Status": "OK"})

data.append({"Parâmetro": "Número de tripulantes",
             "Requisito": "Definido por peso",
             "Valor": n_tripulantes,
             "Status": "Info"})

data.append({"Parâmetro": "MTOW com todos os PAX",
             "Requisito": "Minimizar",
             "Valor": f"{W0:.2f} N",
             "Status": "Check"})

data.append({"Parâmetro": "Saídas de emergência",
             "Requisito": "Regulamento",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Altura do compartimento de carga",
             "Requisito": "≥ 70 cm",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Posição de início de estol (2y/b)",
             "Requisito": "≤ 2y/b aileron",
             "Valor": None,
             "Status": None})

# Distância de decolagem calculada
data.append({"Parâmetro": "Pista de decolagem (ISA, SL)",
             "Requisito": f"≤ {distance_takeoff_req:.0f} m",
             "Valor": f"{distance_takeoff_calc:.2f} m",
             "Status": "OK" if distance_takeoff_calc <= distance_takeoff_req else "NOK"})

# Distância de pouso calculada
data.append({"Parâmetro": "Pista de pouso (ISA, SL)",
             "Requisito": f"≤ {distance_landing_req:.0f} m",
             "Valor": f"{distance_landing_calc:.2f} m",
             "Status": "OK" if distance_landing_calc <= distance_landing_req else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.111)",
             "Requisito": "≥ 1.2%",
             "Valor": f"{grad_25_111:.2f} %",
             "Status": "OK" if grad_25_111 >= 1.2 else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - gear down)",
             "Requisito": "≥ 0.0%",
             "Valor": f"{grad_121_a:.2f} %",
             "Status": "OK" if grad_121_a >= 0.0 else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - gear up)",
             "Requisito": "≥ 2.4%",
             "Valor": f"{grad_121_b:.2f} %",
             "Status": "OK" if grad_121_b >= 2.4 else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - flaps up)",
             "Requisito": "≥ 1.2%",
             "Valor": f"{grad_121_c:.2f} %",
             "Status": "OK" if grad_121_c >= 1.2 else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.119)",
             "Requisito": "≥ 3.2%",
             "Valor": f"{grad_25_119:.2f} %",
             "Status": "OK" if grad_25_119 >= 3.2 else "NOK"})

data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - arremetida)",
             "Requisito": "≥ 2.1%",
             "Valor": f"{grad_121_d:.2f} %",
             "Status": "OK" if grad_121_d >= 2.1 else "NOK"})

data.append({"Parâmetro": "Mach de cruzeiro",
             "Requisito": "0.75 - 0.80",
             "Valor": f"{Mach_cruise:.2f}",
             "Status": "OK"})

data.append({"Parâmetro": "Altitude de cruzeiro",
             "Requisito": "≥ 35000 ft",
             "Valor": f"{altitude_cruise/0.3048:.2f} ft",
             "Status": "OK" if altitude_cruise/0.3048 >= 35000 else "NOK"})

data.append({"Parâmetro": "Excesso de tração na condição de cruzeiro",
             "Requisito": "≥ 0.0",
             "Valor": f"{excess_thrust_cruise:.2f} N (Ref T0 req)",
             "Status": "OK" if excess_thrust_cruise >= 0 else "NOK"})

data.append({"Parâmetro": "Margem estática dianteira",
             "Requisito": "≤ 30%",
             "Valor": f"{SM_fwd*100:.2f} %",
             "Status": "OK" if SM_fwd*100 <= 30 else "NOK"})

data.append({"Parâmetro": "Margem estática traseira",
             "Requisito": "≥ 5%",
             "Valor": f"{SM_aft*100:.2f} %",
             "Status": "OK" if SM_aft*100 >= 5 else "NOK"})

data.append({"Parâmetro": "Estabilidade estática direcional (Cn_beta)",
             "Requisito": "≥ 0",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade estática lateral (Cl_beta)",
             "Requisito": "≤ 0",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade dinâmica - período curto",
             "Requisito": "≤ Nível II",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade dinâmica - fugoidal",
             "Requisito": "≤ Nível II",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade dinâmica - espiral",
             "Requisito": "≤ Nível II",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade dinâmica - dutch roll",
             "Requisito": "≤ Nível II",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Estabilidade dinâmica - rolamento",
             "Requisito": "≤ Nível II",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Deflexão máxima de profundor",
             "Requisito": "≤ 30º",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "CL da empenagem vertical para condição OEI",
             "Requisito": "≤ 0.75",
             "Valor": None,
             "Status": None})

data.append({"Parâmetro": "Fração máxima de peso na bequilha",
             "Requisito": "≤ 18%",
             "Valor": f"{max(frac_nlg_fwd, frac_nlg_aft)*100:.2f} %",
             "Status": "OK" if max(frac_nlg_fwd, frac_nlg_aft)*100 <= 18 else "NOK"})

data.append({"Parâmetro": "Fração mínima de peso na bequilha",
             "Requisito": "≥ 5%",
             "Valor": f"{min(frac_nlg_fwd, frac_nlg_aft)*100:.2f} %",
             "Status": "OK" if min(frac_nlg_fwd, frac_nlg_aft)*100 >= 5 else "NOK"})

data.append({"Parâmetro": "Ângulo de tipback",
             "Requisito": "≥ 15º",
             "Valor": f"{np.degrees(alpha_tipback):.2f} º",
             "Status": "OK" if np.degrees(alpha_tipback) >= 15 else "NOK"})

data.append({"Parâmetro": "Ângulo de tailstrike",
             "Requisito": "≥ 10º",
             "Valor": f"{np.degrees(alpha_tailstrike):.2f} º",
             "Status": "OK" if np.degrees(alpha_tailstrike) >= 10 else "NOK"})

data.append({"Parâmetro": "Ângulo de overturn",
             "Requisito": "≤ 63º",
             "Valor": f"{np.degrees(phi_overturn):.2f} º",
             "Status": "OK" if np.degrees(phi_overturn) <= 63 else "NOK"})

data.append({"Parâmetro": "Fração da envergadura ocupada pelo tanque",
             "Requisito": "≤ 95%",
             "Valor": f"{b_tank_b_w*100:.2f} %",
             "Status": "OK" if b_tank_b_w*100 <= 95 else "NOK"})

data.append({"Parâmetro": "Envergadura máxima (Gate Code C)",
             "Requisito": "≤ 36 m",
             "Valor": f"{b_wing:.2f} m",
             "Status": "OK" if b_wing <= 36 else "NOK"})

data.append({"Parâmetro": "Motores livres do spray gerado por TDP",
             "Requisito": "Sim",
             "Valor": None,
             "Status": None})

df = pd.DataFrame(data)

output_file = Path(__file__).parent.parent / 'Resultados' / 'Compliance_Matrix.xlsx'
output_file.parent.mkdir(parents=True, exist_ok=True)
df.to_excel(output_file, index=False)

print(df[['Parâmetro', 'Valor', 'Status']].to_string())

print("\nDetalhes da Análise de Peso e Balanceamento:\n")
print(f"W0 = {W0}")
print(f"Wf = {Wf}\n")
print(f"T0 = {T0}\n")
print(f"deltaS_wlan = {deltaS_wlan}\n")
print(f"SM_fwd = {SM_fwd}")
print(f"SM_aft = {SM_aft}\n")
print(f"b_tank_b_w = {b_tank_b_w}\n")
print(f"frac_nlg_fwd = {frac_nlg_fwd}")
print(f"frac_nlg_aft = {frac_nlg_aft}")
print(f"alpha_tipback = {alpha_tipback}")
print(f"alpha_tailstrike = {alpha_tailstrike}")
print(f"phi_overturn = {phi_overturn}")

print(f"\nDistância de decolagem calculada (ISA, SL): {distance_takeoff_calc:.2f} m")
print(f"Distância de pouso calculada (ISA, SL): {distance_landing_calc:.2f} m")

plot3d(aircraft)
