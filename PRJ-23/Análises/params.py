import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import aux_tools as at
import design_tools_optmized as dt
aircraft = dt.my_aircraft()

gravity = 9.81
g = 9.81 # Aceleração da gravidade
T0_guess = 125600 #Chute inicial
W0_guess = 490000.0 #Chute inicial
altitude_cruise = 11000.0000 #Carteado
Mach_cruise = 0.7500000 #Range de 0.75 a 0.80
range_cruise = 3700e3 # Req projeto
range_altcruise = 370400 # 200 NM
loiter_time = 2700.00000 # 45 minutos
altitude_altcruise = 4572.00000 # Caso de testess
Mach_altcruise = 0.40000000 #Carteado

distance_takeoff = 1800.0 #Req projeto
distance_landing = 1150.0 # Req projeto

TO_flap_def = 20 * np.pi / 180 
LD_flap_def = aircraft['data']['flap']['max_def']
TO_slat_def = 0
LD_slat_def = 0
h_ground = 10.668
altitude_cruise = 11000
altitude_takeoff = 0.0
altitude_landing = 0.0
MLW_frac = 0.84

distance_takeoff = 1800.0 #Req projeto
distance_landing = 1150.0 # Req projeto


new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave

W0, Wf, T0, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

Ixx,Iyy,Izz,Ixy,Ixz,Iyz = dt.moment_of_inertia(aircraft)

print(f"Ixx = {Ixx}")
print(f"Iyy = {Iyy}")
print(f"Izz = {Izz}")
print(f"Ixy = {Ixy}")
print(f"Ixz = {Ixz}")
print(f"Iyz = {Iyz}")
print()

T, p, rho, mi = at.atmosphere(z=altitude_cruise, Tba=288.15)
a = np.sqrt(1.4 * 287.05 * T)  # Speed of sound [m/s]

# Create parameter table
table_data = {
    'Parâmetro': ['S_ref', 'c_ref', 'b_ref', 'm', 'I_xx', 'I_yy', 'I_zz', 'I_xz', 'i_p', 'x_p', 'z_p', 'T_max', 'V', 'h'],
    'Explicação': [
        'área de referência [m]',
        'corda de referência [m]',
        'envergadura de referência [m]',
        'massa da aeronave [kg]',
        'momento de inércia [kg·m²]',
        'momento de inércia [kg·m²]',
        'momento de inércia [kg·m²]',
        'momento de inércia [kg·m²]',
        'incidência do motor [°]',
        'posição longitudinal do motor [m]',
        'posição vertical do motor [m]',
        'tração máxima [N]',
        'velocidade de voo [m/s]',
        'altitude de voo [m]'
    ],
    'Valor': [
        aircraft['geo_param']['wing']['S'],
        aircraft['dimensions']['wing']['cm'],
        aircraft['dimensions']['wing']['b'],
        W0 / gravity,
        f"{Ixx:.2f}",
        f"{Iyy:.2f}",
        f"{Izz:.2f}",
        f"{Ixz:.2f}",
        0,
        aircraft['dimensions']['nacelle']['xn'],
        aircraft['dimensions']['nacelle']['zn'],
        # 116049.63578638842,
        None,
        Mach_cruise * a,
        altitude_cruise
    ],
    'Fonte': ['designTool'] * 14
}

df = pd.DataFrame(table_data)

aero, CLmax = dt.aerodynamics(aircraft, Mach_cruise, altitude_cruise, 0, LD_flap_def, LD_slat_def, 0, 0, W0_guess)

print(aero)

print(df)

df.to_excel('PRJ-23\\Resultados\\aircraft_parameters.xlsx', index=False)

