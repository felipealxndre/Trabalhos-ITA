#Grupo - Pinguins de Madagascar
#Integrante:
# - Thiago Akira Missato

# Importar os arquivos
import design_tools_default as dt
import aux_tools as at
import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt
from scipy.stats import linregress

aircraft = dt.my_aircraft() #Define a aeronave nossa.
new_dimensions = dt.geometry(aircraft) #Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave.

#Caso de testes



gravity = 9.81
W0_guess = 43090 * gravity
T0_guess = 125600

Mach_cruise = 0.77
altitude_cruise = 11000
range_cruise = 2390000.0

Mach_altcruise = 0.4
range_altcruise = 370000
altitude_altcruise = 4572

loiter_time = 2700

altitude_takeoff = 0
distance_takeoff = 1520
TO_flap_def = 0.34906585039887
TO_slat_def = 0

altitude_landing = 0
distance_landing = 1520
LD_flap_def = 0.69813170079773
LD_slat_def = 0
MLW_frac = 0.84

# referencias

W0_ref, Wf_ref, We_ref, deltaS_wlan, SM_fwd_ref, SM_aft_ref, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise,Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

references = {
    'W0': W0_ref,
    'Wf': Wf_ref,
    'We': We_ref,
    'SM_fwd': SM_fwd_ref,
    'SM_aft': SM_aft_ref,
    'deltaS_wlan': deltaS_wlan,
    'b_tank_b_w': b_tank_b_w,
    'frac_nlg_fwd': frac_nlg_fwd,
    'frac_nlg_aft': frac_nlg_aft,
    'alpha_tipback': alpha_tipback,
    'alpha_tailstrike': alpha_tailstrike,
    'phi_overturn': phi_overturn
}


# Definição de parâmetros de input

# asa
Sw = aircraft['geo_param']['wing']['S']
AR_w = aircraft['geo_param']['wing']['AR']
sweep_w = aircraft['geo_param']['wing']['sweep']
delta_w = aircraft['geo_param']['wing']['dihedral']
xr_w = aircraft['geo_param']['wing']['xr']
tcr_w = aircraft['geo_param']['wing']['tcr']
tct_w = aircraft['geo_param']['wing']['tct']

# empenagens
Cht = aircraft['geo_param']['EH']['Cht']
Lc = aircraft['geo_param']['EH']['Lc']
Cvt = aircraft['geo_param']['EV']['Cvt']
Lb = aircraft['geo_param']['EV']['Lb']

# misc
x_tailstrike = aircraft['data']['misc']['x_tailstrike']
z_tailstrike = aircraft['data']['misc']['z_tailstrike']


# landing gear
xnlg = aircraft['dimensions']['ldg']['xnlg']
xmlg = aircraft['dimensions']['ldg']['xmlg']
ymlg = aircraft['dimensions']['ldg']['ymlg']


input_params = [Sw, AR_w, sweep_w, delta_w, xr_w, tcr_w, tct_w, Cht, Lc, Cvt, Lb, x_tailstrike, z_tailstrike, xnlg, xmlg, ymlg, Mach_cruise, range_cruise]
param_names = ['Sw', 'AR_w', 'sweep_w', 'delta_w', 'xr_w', 'tcr_w', 'tct_w', 'Cht', 'Lc', 'Cvt', 'Lb', 'x_tailstrike', 'z_tailstrike', 'xnlg', 'xmlg', 'ymlg', 'Mach_cruise', 'range_cruise']
output = []
var_index = 0

for i, param in enumerate(input_params):
   
    aircraft = dt.my_aircraft()
    current_mach = 0.77
    current_range = 2390000.0
    
    if i == 0:  # Sw
        aircraft['geo_param']['wing']['S'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 1:  # AR_w
        aircraft['geo_param']['wing']['AR'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 2:  # sweep_w
        aircraft['geo_param']['wing']['sweep'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 3:  # delta_w
        aircraft['geo_param']['wing']['dihedral'] = param + 2 * np.pi / 180
        sensitivity_param = (2 * np.pi / 180) / param
    elif i == 4:  # xr_w
        aircraft['geo_param']['wing']['xr'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 5:  # tcr_w
        aircraft['geo_param']['wing']['tcr'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 6:  # tct_w
        aircraft['geo_param']['wing']['tct'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 7:  # Cht
        aircraft['geo_param']['EH']['Cht'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 8:  # Lc
        aircraft['geo_param']['EH']['Lc'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 9:  # Cvt
        aircraft['geo_param']['EV']['Cvt'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 10:  # Lb
        aircraft['geo_param']['EV']['Lb'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 11:  # x_tailstrike
        aircraft['data']['misc']['x_tailstrike'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 12:  # z_tailstrike
        aircraft['data']['misc']['z_tailstrike'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 13:  # xnlg
        aircraft['dimensions']['ldg']['xnlg'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 14:  # xmlg
        aircraft['dimensions']['ldg']['xmlg'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 15:  # ymlg
        aircraft['dimensions']['ldg']['ymlg'] = param * 1.02
        sensitivity_param = 0.02
    elif i == 16:  # Mach_cruise
        current_mach = param + 0.02
        sensitivity_param = 0.02 / param
    elif i == 17:  # range_cruise
        current_range = param * 1.02
        sensitivity_param = 0.02
    else:
        current_range = param
        sensitivity_param = 0.02 / param

    new_dimensions = dt.geometry(aircraft) #Calcula as novas dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave.
    W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, current_mach, altitude_cruise, current_range,Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
    
    sensibilities = {
        'W0': (W0 - references['W0']) / references['W0'],
        'Wf': (Wf - references['Wf']) / references['Wf'],
        'We': (We - references['We']) / references['We'],
        'SM_fwd': (SM_fwd - references['SM_fwd']) / references['SM_fwd'],
        'SM_aft': (SM_aft - references['SM_aft']) / references['SM_aft'],
        'deltaS_wlan': (deltaS_wlan - references['deltaS_wlan']) / references['deltaS_wlan'],
        'b_tank_b_w': (b_tank_b_w - references['b_tank_b_w']) / references['b_tank_b_w'],
        'frac_nlg_fwd': (frac_nlg_fwd - references['frac_nlg_fwd']) / references['frac_nlg_fwd'],
        'frac_nlg_aft': (frac_nlg_aft - references['frac_nlg_aft']) / references['frac_nlg_aft'],
        'alpha_tipback': (alpha_tipback - references['alpha_tipback']) / references['alpha_tipback'],
        'alpha_tailstrike': (alpha_tailstrike - references['alpha_tailstrike']) / references['alpha_tailstrike'],
        'phi_overturn': (phi_overturn - references['phi_overturn']) / references['phi_overturn']
    }

    current_output = []
    for key in sensibilities.keys():
        current_output.append(round(sensibilities[key]/sensitivity_param, 2))

    print(f"Parameter: {param_names[i]}, Sensitivities: {current_output}")
    output.append(current_output)

# Convert output to DataFrame
sensitivity_columns = ['W0', 'Wf', 'We', 'SM_fwd', 'SM_aft', 'deltaS_wlan', 'b_tank_b_w', 
                      'frac_nlg_fwd', 'frac_nlg_aft', 'alpha_tipback', 'alpha_tailstrike', 'phi_overturn']

sensibilidade = pd.DataFrame(output, index=param_names, columns=sensitivity_columns)

sensibilidade.to_excel('sensibilidade.xlsx')






