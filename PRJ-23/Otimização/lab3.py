# IMPORTS
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from pprint import pprint
import copy  # Add this import
# removing warnings in terminal
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import design_tools_optmized as dt
import aux_tools as at
import auxmod as am

# input fixed parameters
gravity = 9.81
W0_guess = 43090 * gravity
T0_guess = 125600

Mach_cruise = 0.75
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

aircraft = dt.my_aircraft()

# updating from the previously optimization - Point B 

# # Updated aircraft parameters - review
# aircraft['geo_param']['wing']['AR'] = 11.97
# aircraft['geo_param']['wing']['S'] = 98.77
# aircraft['geo_param']['wing']['sweep'] = 0.23
# aircraft['geo_param']['EH']['Cht'] = 0.80
# aircraft['dimensions']['ldg']['xmlg'] = 15.42
# aircraft['dimensions']['ldg']['ymlg'] = 4.11
# aircraft['dimensions']['ldg']['xnlg'] = 1.0



new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave




#aero, CLmax = dt.aerodynamics(aircraft=aircraft, Mach=Mach_cruise, altitude=altitude_cruise, n_engines_failed=0, flap_def=LD_flap_def, slat_def=LD_slat_def, lg_down=0, h_ground=0, W0_guess=W0_guess)



W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)



# --- Create DataFrame for reference section data ---

# Atmospheric properties at cruise altitude
T, p, rho, mi = at.atmosphere(z=altitude_cruise, Tba=288.15)
a = np.sqrt(1.4 * 287.05 * T)  # Speed of sound [m/s]

# Reference parameters
h = altitude_cruise
W = W0  # Use the optimized weight
S_ref = aircraft['geo_param']['wing']['S']
a_inf = a
c_ref = aircraft['dimensions']['wing']['cm']
cl_ref = 0.4342478517869541  # CL at project point
Re_ref = (rho * Mach_cruise * a_inf * c_ref) / mi

#tc_ref = aircraft['geo_param']['wing']['tc']

print('tcr:', aircraft['geo_param']['wing']['tcr'])
print('tct:', aircraft['geo_param']['wing']['tct'])
tc_ref = (aircraft['geo_param']['wing']['tcr'] + aircraft['geo_param']['wing']['tct'])*0.5

# Build DataFrame
df_ref = pd.DataFrame({
    'Parâmetro': [
        'h', 'M_∞', 'W', 'S_ref', 'rho_∞', 'a_∞', 'c_ref', 'c_{l_ref}', 'Re_{ref}', '(t/c)_{ref}'
    ],
    'Valor': [
        h, Mach_cruise, W, S_ref, rho, a_inf, c_ref, cl_ref, Re_ref, tc_ref
    ],
    'Descrição': [
        'Altitude do ponto de projeto [ft]',
        'Mach do ponto de projeto',
        'Peso da aeronave no ponto de projeto [kgf]',
        'Área de referência da aeronave [m²]',
        'Densidade do ar no ponto de projeto [kg/m³]',
        'Velocidade do som no ponto de projeto [m/s]',
        'Corda da seção de referência [m]',
        'Coeficiente de sustentação da seção de referência',
        'Reynolds para a seção de referência',
        'Espessura relativa base da seção de referência'
    ]
})




print(df_ref)

df_ref.to_excel('Otimização\\Resultados\\df_ref.xlsx', index=False)


