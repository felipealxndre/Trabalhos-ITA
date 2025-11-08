import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import design_tools as dt
from aux_tools import atmosphere
import seaborn as sns


gravity = 9.81
W0_guess = 43090*gravity
T0_guess = 125600
Mach_cruise = 0.75
altitude_cruise = 11000
range_cruise  = 2390000
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

# configuração inicial da aeronave
aircraft = dt.my_aircraft()

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

W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)


Ixx,Iyy,Izz,Ixy,Ixz,Iyz = dt.moment_of_inertia(aircraft)
# saving in txt file
with open(os.path.join('aircraft_data.txt'), 'w') as f:
    f.write(f"""

aircraft.S = {aircraft['geo_param']['wing']['S']}; % Reference wing area [m^2]
aircraft.c = {aircraft['dimensions']['wing']['cm']}; % Wing mean aerodynamic chord [m]
aircraft.b = {aircraft['dimensions']['wing']['b']}; % Wing span [m]
%--------------------------------------------------------------------------
% Inertia properties:
aircraft.m = {W0/9.81}; % Mass [kg]
aircraft.Ixx = {Ixx}; % I_xx [kg.m^2]
aircraft.Iyy = {Iyy}; % I_yy [kg.m^2]
aircraft.Izz = {Izz}; % I_zz [kg.m^2]
aircraft.Ixz = {Ixz}; % I_xz [kg.m^2]


%--------------------------------------------------------------------------
% pegar do enunciado
% Propulsive model properties:
aircraft.iota_1_deg = 2.0; % Left engine incidence [deg] 
aircraft.tau_1_deg = 1.5; % Left engine toe-in angle [deg] 
aircraft.x_1 = 12.0; % Left engine thrust action point x_b coordinate [m] (de design_tools.py, 'xn')
aircraft.y_1 = -4.2; % Left engine thrust action point y_b coordinate [m] (de design_tools.py, -'yn')
aircraft.z_1 = -2.5192; % Left engine thrust action point z_b coordinate [m] (de design_tools.py, 'zn')
aircraft.iota_2_deg = 2.0; % Right engine incidence [deg] 
aircraft.tau_2_deg = -1.5; % Right engine toe-in angle [deg] 
aircraft.x_2 = 12.0; % Right engine thrust action point x_b coordinate [m] (de design_tools.py, 'xn')
aircraft.y_2 = 4.2; % Right engine thrust action point y_b coordinate [m] (de design_tools.py, 'yn')
aircraft.z_2 = -2.5192; % Right engine thrust action point z_b coordinate [m] (de design_tools.py, 'zn')
aircraft.Tmax = 100000; % Maximum sea-level thrust for a single engine [N] 
aircraft.n_rho = 0.8; % Air density exponent in the thrust model [-]

%--------------------------------------------------------------------------
% outros parametros

trim_par = struct('V',{Mach_cruise*np.sqrt(1.4 * 8.31 * atmosphere(z= altitude_cruise, Tba=288.15)[0])},'h',{altitude_cruise},'gamma_deg',{-3},...
    'thetadot_deg_s',0);


% dados lab 4

Cdalpha = 0.002419;
Cdalha2 = 0.031170;
""")
    