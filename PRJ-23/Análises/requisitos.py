import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import design_tools as dt
import aux_tools as at

aircraft = dt.my_aircraft()

# Updated aircraft parameters - review
aircraft['geo_param']['wing']['AR'] = 9.151089
aircraft['geo_param']['wing']['S'] = 99.033744
aircraft['geo_param']['wing']['sweep'] = 0.230809
aircraft['geo_param']['EH']['Cht'] = 0.900003
aircraft['dimensions']['ldg']['xnlg'] = 1.0
aircraft['dimensions']['ldg']['xmlg'] = 15.610284
aircraft['dimensions']['ldg']['ymlg'] = 2.470000 


new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)

g = 9.81
# T0_guess = 125600
# W0_guess = 490000.0

# range_cruise = 3700e3
# distance_takeoff = 1800.0
# distance_landing = 1150.0
# Mach_cruise = 0.77
# altitude_cruise = 11000.0
# loiter_time = 45 * 60
# range_altcruise = 200 * 1852
# altitude_altcruise = 4572.0 
# Mach_altcruise = 0.40
# altitude_takeoff = 0.0
# altitude_landing = 0.0



# TO_flap_def = 20 * np.pi / 180
# TO_slat_def = 0
# LD_flap_def = aircraft['data']['flap']['max_def']
# LD_slat_def = 0
# h_ground = 10.668
# MLW_frac = 0.84

W0_guess = 422712.90000000002328
T0_guess = 125600
TO_flap_def = 0.34906585039887
LD_flap_def = 0.69813170079773
TO_slat_def = 0.0
LD_slat_def = 0.0
h_ground = 10.668
altitude_cruise = 11000
Mach_cruise = 0.77
range_cruise = 3700e3

loiter_time = 45*60
altitude_altcruise = 4572
Mach_altcruise = 0.4
range_altcruise = 200*1852
altitude_takeoff = 0.0
distance_takeoff = 1800
altitude_landing = 0.0
distance_landing = 1150
MLW_frac = 0.84


W0, Wf, T0, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, \
frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = \
    dt.analyze(aircraft, W0_guess, T0_guess, 
                Mach_cruise, altitude_cruise, range_cruise,
                Mach_altcruise, range_altcruise, altitude_altcruise, 
                loiter_time, 
                altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, 
                altitude_landing, distance_landing, LD_flap_def, LD_slat_def, 
                MLW_frac)

aircraft['weights']['W0'] = W0
T0_available, T0vec_req, S_wlan_req = dt.performance(aircraft, 
                                    TO_flap_def, LD_flap_def, 
                                    TO_slat_def, LD_slat_def,
                                    h_ground, 
                                    altitude_takeoff, distance_takeoff, 
                                    altitude_landing, distance_landing, MLW_frac, 
                                    altitude_cruise, Mach_cruise)

b_wing = aircraft['dimensions']['wing']['b']
n_pax = 90
n_tripulantes = 5
excess_thrust = T0_available - T0vec_req
excess_grad = excess_thrust / W0

grad_25_111 = 1.2 + (excess_grad[2] * 100)
grad_121_a = 0.0 + (excess_grad[3] * 100)
grad_121_b = 2.4 + (excess_grad[4] * 100)
grad_121_c = 1.2 + (excess_grad[5] * 100)
grad_25_119 = 3.2 + ((excess_thrust[6] / (W0 * MLW_frac)) * 100)
grad_121_d = 2.1 + ((excess_thrust[7] / (W0 * MLW_frac)) * 100)

excess_thrust_cruise = T0_available - T0vec_req[1]


data = []

data.append({"Parâmetro": "Número de passageiros", "Requisito": "90", "Valor": n_pax, "Status": "OK"})
data.append({"Parâmetro": "Número de tripulantes", "Requisito": "Definido por peso", "Valor": n_tripulantes, "Status": "Info"})
data.append({"Parâmetro": "MTOW com todos os PAX", "Requisito": "Minimizar", "Valor": f"{W0:.2f} N", "Status": "Check"})
data.append({"Parâmetro": "Saídas de emergência", "Requisito": "Regulamento", "Valor": None, "Status": None})
data.append({"Parâmetro": "Altura do compartimento de carga", "Requisito": "≥ 70 cm", "Valor": None, "Status": None})
data.append({"Parâmetro": "Posição de início de estol (2y/b)", "Requisito": "≤ 2y/b aileron", "Valor": None, "Status": None})
data.append({"Parâmetro": "Pista de decolagem (ISA, SL)", "Requisito": f"≤ {distance_takeoff:.0f} m", "Valor": f"{distance_takeoff:.2f} m", "Status": "OK"})
data.append({"Parâmetro": "Pista de pouso (ISA, SL)", "Requisito": f"≤ {distance_landing:.0f} m", "Valor": f"{distance_landing:.2f} m", "Status": "OK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.111)", "Requisito": "≥ 1.2%", "Valor": f"{grad_25_111:.2f} %", "Status": "OK" if grad_25_111 >= 1.2 else "NOK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - gear down)", "Requisito": "≥ 0.0%", "Valor": f"{grad_121_a:.2f} %", "Status": "OK" if grad_121_a >= 0.0 else "NOK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - gear up)", "Requisito": "≥ 2.4%", "Valor": f"{grad_121_b:.2f} %", "Status": "OK" if grad_121_b >= 2.4 else "NOK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - flaps up)", "Requisito": "≥ 1.2%", "Valor": f"{grad_121_c:.2f} %", "Status": "OK" if grad_121_c >= 1.2 else "NOK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.119)", "Requisito": "≥ 3.2%", "Valor": f"{grad_25_119:.2f} %", "Status": "OK" if grad_25_119 >= 3.2 else "NOK"})
data.append({"Parâmetro": "Gradiente de subida (FAR 25.121 - arremetida)", "Requisito": "≥ 2.1%", "Valor": f"{grad_121_d:.2f} %", "Status": "OK" if grad_121_d >= 2.1 else "NOK"})
data.append({"Parâmetro": "Mach de cruzeiro", "Requisito": "0.75 - 0.80", "Valor": f"{Mach_cruise:.2f}", "Status": "OK"})
data.append({"Parâmetro": "Altitude de cruzeiro", "Requisito": "≥ 35000 ft", "Valor": f"{altitude_cruise/0.3048:.2f} ft", "Status": "OK" if altitude_cruise/0.3048 >= 35000 else "NOK"})
data.append({"Parâmetro": "Excesso de tração na condição de cruzeiro", "Requisito": "≥ 0.0", "Valor": f"{excess_thrust_cruise:.2f} N (Ref T0 req)", "Status": "OK" if excess_thrust_cruise >= 0 else "NOK"})
data.append({"Parâmetro": "Margem estática dianteira", "Requisito": "≤ 30%", "Valor": f"{SM_fwd*100:.2f} %", "Status": "OK" if SM_fwd*100 <= 30 else "NOK"})
data.append({"Parâmetro": "Margem estática traseira", "Requisito": "≥ 5%", "Valor": f"{SM_aft*100:.2f} %", "Status": "OK" if SM_aft*100 >= 5 else "NOK"})
data.append({"Parâmetro": "Estabilidade estática direcional (Cn_beta)", "Requisito": "≥ 0", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade estática lateral (Cl_beta)", "Requisito": "≤ 0", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade dinâmica - período curto", "Requisito": "≤ Nível II", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade dinâmica - fugoidal", "Requisito": "≤ Nível II", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade dinâmica - espiral", "Requisito": "≤ Nível II", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade dinâmica - dutch roll", "Requisito": "≤ Nível II", "Valor": None, "Status": None})
data.append({"Parâmetro": "Estabilidade dinâmica - rolamento", "Requisito": "≤ Nível II", "Valor": None, "Status": None})
data.append({"Parâmetro": "Deflexão máxima de profundor", "Requisito": "≤ 30º", "Valor": None, "Status": None})
data.append({"Parâmetro": "CL da empenagem vertical para condição OEI", "Requisito": "≤ 0.75", "Valor": None, "Status": None})
data.append({"Parâmetro": "Fração máxima de peso na bequilha", "Requisito": "≤ 18%", "Valor": f"{max(frac_nlg_fwd, frac_nlg_aft)*100:.2f} %", "Status": "OK" if max(frac_nlg_fwd, frac_nlg_aft)*100 <= 18 else "NOK"})
data.append({"Parâmetro": "Fração mínima de peso na bequilha", "Requisito": "≥ 5%", "Valor": f"{min(frac_nlg_fwd, frac_nlg_aft)*100:.2f} %", "Status": "OK" if min(frac_nlg_fwd, frac_nlg_aft)*100 >= 5 else "NOK"})
data.append({"Parâmetro": "Ângulo de tipback", "Requisito": "≥ 15º", "Valor": f"{np.degrees(alpha_tipback):.2f} º", "Status": "OK" if np.degrees(alpha_tipback) >= 15 else "NOK"})
data.append({"Parâmetro": "Ângulo de tailstrike", "Requisito": "≥ 10º", "Valor": f"{np.degrees(alpha_tailstrike):.2f} º", "Status": "OK" if np.degrees(alpha_tailstrike) >= 10 else "NOK"})
data.append({"Parâmetro": "Ângulo de overturn", "Requisito": "≤ 63º", "Valor": f"{np.degrees(phi_overturn):.2f} º", "Status": "OK" if np.degrees(phi_overturn) <= 63 else "NOK"})
data.append({"Parâmetro": "Fração da envergadura ocupada pelo tanque", "Requisito": "≤ 95%", "Valor": f"{b_tank_b_w*100:.2f} %", "Status": "OK" if b_tank_b_w*100 <= 95 else "NOK"})
data.append({"Parâmetro": "Envergadura máxima (Gate Code C)", "Requisito": "≤ 36 m", "Valor": f"{b_wing:.2f} m", "Status": "OK" if b_wing <= 36 else "NOK"})
data.append({"Parâmetro": "Motores livres do spray gerado por TDP", "Requisito": "Sim", "Valor": None, "Status": None})

df = pd.DataFrame(data)

output_file = Path(__file__).parent.parent / 'Resultados' / 'Compliance_Matrix.xlsx'
df.to_excel(output_file, index=False)

print("-" * 83)
print(df[['Parâmetro', 'Valor', 'Status']].to_string())


print("\nDetalhes do Análise de Peso e Balanceamento:\n")
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