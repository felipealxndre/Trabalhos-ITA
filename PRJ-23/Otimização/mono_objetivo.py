# IMPORTS
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import seaborn as sns  # Add seaborn import

from pprint import pprint
from time import time
import copy  # Add this import
# removing warnings in terminal
import warnings
import pandas as pd

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import auxmod as am
import design_tools as dt
import aux_tools as at
from plot3d import plot3d

warnings.filterwarnings("ignore")


# EXECUTION

# Define history list
Xlist = []  # design variables points
flist = []  # objective function values
glist = []  # constraints values


# Define objective function
def objfun(X):

    # design variables - Sw, ARw, sweep, Cht, X_nlg, X_mlg, Y_mlg

    # desnormalizando
    AR_w = X[0] * (bounds['AR_w'][1] - bounds['AR_w'][0]) + bounds['AR_w'][0]
    S_w = X[1] * (bounds['Sw'][1] - bounds['Sw'][0]) + bounds['Sw'][0]
    sweep_w = (X[2] * (bounds['sweep_w'][1] - bounds['sweep_w'][0]) + bounds['sweep_w'][0]) *  np.pi/180
    Cht = X[3] * (bounds['Cht'][1] - bounds['Cht'][0]) + bounds['Cht'][0]
    xnlg = X[4] * (bounds['xnlg'][1] - bounds['xnlg'][0]) + bounds['xnlg'][0]
    xmlg = X[5] * (bounds['xmlg'][1] - bounds['xmlg'][0]) + bounds['xmlg'][0]
    ymlg = X[6] * (bounds['ymlg'][1] - bounds['ymlg'][0]) + bounds['ymlg'][0]

    # inputing other restrictions
    aircraft['geo_param']['wing']['AR'] = AR_w
    aircraft['geo_param']['wing']['S'] = S_w
    aircraft['geo_param']['wing']['sweep'] = sweep_w
    aircraft['geo_param']['EH']['Cht'] = Cht
    aircraft['dimensions']['ldg']['xnlg'] = xnlg
    aircraft['dimensions']['ldg']['xmlg'] = xmlg
    aircraft['dimensions']['ldg']['ymlg'] = ymlg

    new_dimensions = dt.geometry(aircraft) #Calcula as dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave
    # Compute objective function

    print('-'*50)
    print('History:')
    print('AR:', AR_w)
    print('Sw:', S_w)
    print('Sweep:', sweep_w)
    print('Cht:', Cht)
    print('xnlg:', xnlg)
    print('xmlg:', xmlg)
    print('ymlg:', ymlg)

    # analyse to see objective functions
    W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
    
    f = W0 / references['W0']

    # Save results for the plot
    Xlist.append(X)
    flist.append(W0)

    return f

# Define constraints
def confun(X):

    # design variables - Sw, ARw, sweep, Cht, xnlg, xmlg, ymlg

    # desnormalizando
    AR_w = X[0] * (bounds['AR_w'][1] - bounds['AR_w'][0]) + bounds['AR_w'][0]
    S_w = X[1] * (bounds['Sw'][1] - bounds['Sw'][0]) + bounds['Sw'][0]
    sweep_w = (X[2] * (bounds['sweep_w'][1] - bounds['sweep_w'][0]) + bounds['sweep_w'][0]) * np.pi/180
    Cht = X[3] * (bounds['Cht'][1] - bounds['Cht'][0]) + bounds['Cht'][0]
    xnlg = X[4] * (bounds['xnlg'][1] - bounds['xnlg'][0]) + bounds['xnlg'][0]
    xmlg = X[5] * (bounds['xmlg'][1] - bounds['xmlg'][0]) + bounds['xmlg'][0]
    ymlg = X[6] * (bounds['ymlg'][1] - bounds['ymlg'][0]) + bounds['ymlg'][0]

    # inputing other restrictions
    aircraft['geo_param']['wing']['AR'] = AR_w
    aircraft['geo_param']['wing']['S'] = S_w
    aircraft['geo_param']['wing']['sweep'] = sweep_w
    aircraft['geo_param']['EH']['Cht'] = Cht
    aircraft['dimensions']['ldg']['xnlg'] = xnlg
    aircraft['dimensions']['ldg']['xmlg'] = xmlg
    aircraft['dimensions']['ldg']['ymlg'] = ymlg

    print('-'*50)
    print('History:')
    print('AR:', AR_w)
    print('Sw:', S_w)
    print('Sweep:', sweep_w)
    print('Cht:', Cht)
    print('xnlg:', xnlg)
    print('xmlg:', xmlg)
    print('ymlg:', ymlg)

    new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave
    
    # Compute objective function

    # analyse to see objective functions
    W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
    # Constraints
    constraints = [
        deltaS_wlan / references['deltaS_wlan'],                # deltaS_wlan >= 0
        (0.3 - SM_fwd) / 0.3,                                   # SM_fwd <= 0.30
        (SM_aft - 0.05) / 0.05,                                 # SM_aft >= 0.05
        (0.18 - frac_nlg_fwd) / 0.18,                           # frac_nlg_fwd <= 0.18
        (frac_nlg_aft - 0.03) / 0.03,                           # frac_nlg_aft >= 0.03
        (alpha_tipback * 180/np.pi - 15) / 15,                  # alpha_tipback >= 15deg 
        (alpha_tailstrike * 180 / np.pi - 10) / 10,             # alpha_tailstrike >= 10deg
        (63 - phi_overturn * 180/np.pi) / 63,                   # phi_overturn <= 63deg
        (0.95 - b_tank_b_w) / 0.95,                             # b_tank_b_w <= 0.95
        (36 - aircraft['dimensions']['wing']['b']) / 36,        # b_w <= 36m
        (4.5 - aircraft['dimensions']['ldg']['ymlg']) / 4.5,    # y_mlg <= 4.5m
        (13.5 - (aircraft['dimensions']['EV']['zt'] - aircraft['dimensions']['ldg']['z'])) / 13.5,   # zt_v - z_mlg <= 13.5m
        
        
        # new constraints
        (aircraft['dimensions']['EV']['L'] - 13) / 13,          # L_h >= 13
        #(0.5 - abs(aircraft['dimensions']['EV']['xm'] - aircraft['dimensions']['EH']['xm'])) / 0.5   # |x_EV - x_EH| <= 0.5m

    ]
    # Optionally, save results for plotting
    glist.append([constraints])


    return constraints

# Create list of constraints

con1 = {
    'type': 'ineq',
    'fun': confun
}

cons = [con1]

# Define starting point and reference values 
aircraft = dt.my_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave.

original_aircraft = copy.deepcopy(aircraft)  # Use deepcopy instead of copy


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

# creating the references - from the reference group aircraft
reference_keys = ['W0', 'Wf', 'We', 'deltaS_wlan', 'SM_fwd', 'SM_aft', 'b_tank_b_w', 'frac_nlg_fwd', 'frac_nlg_aft', 'alpha_tipback', 'alpha_tailstrike', 'phi_overturn']

reference_values = dt.analyze(
    original_aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise,
    Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time,
    altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def,
    altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac
)
references = dict(zip(reference_keys, reference_values))

print('References: ')
pprint(references)



# Define design of the 7 variable bounds - dictonary of lists - [lower, upper]
bounds = {
    'AR_w': [7.1, 12],
    'Sw': [80, 120],
    'sweep_w': [17, 30], # degrees
    'Cht': [1.3, 1.6],
    'xnlg': [2, 3],
    'xmlg': [15, 20],
    'ymlg': [2, 6]
}

# normalizados -> value = (value - lower) /(upper - lower)
X0 = {
    'AR_w': (aircraft['geo_param']['wing']['AR'] - bounds['AR_w'][0])/(bounds['AR_w'][1] - bounds['AR_w'][0]),
    'Sw': (aircraft['geo_param']['wing']['S'] - bounds['Sw'][0])/(bounds['Sw'][1] - bounds['Sw'][0]),
    'sweep_w': (aircraft['geo_param']['wing']['sweep']*180/np.pi - bounds['sweep_w'][0])/(bounds['sweep_w'][1] - bounds['sweep_w'][0]),
    'Cht': (aircraft['geo_param']['EH']['Cht'] - bounds['Cht'][0])/(bounds['Cht'][1] - bounds['Cht'][0]),
    'xnlg': (aircraft['dimensions']['ldg']['xnlg'] - bounds['xnlg'][0])/(bounds['xnlg'][1] - bounds['xnlg'][0]),
    'xmlg': (aircraft['dimensions']['ldg']['xmlg'] - bounds['xmlg'][0])/(bounds['xmlg'][1] - bounds['xmlg'][0]),
    'ymlg': (aircraft['dimensions']['ldg']['ymlg'] - bounds['ymlg'][0])/(bounds['ymlg'][1] - bounds['ymlg'][0]),
}

X0 = np.array(list(X0.values()))
print('Starting Point: \n')
print(X0)



# Additional options
options = {'maxiter':200, 'ftol':1e-6}

start_time = time()

# Run optimizer
result = minimize(objfun, X0,
                  constraints=cons, 
                  bounds= [[0, 1]] * len(bounds), # normalized bounds - all between 0 ad 1
                  method='slsqp', options=options, )

end_time = time()

# Print results
print(result)
print(f"\nTime of optimization: {end_time - start_time}")
Xopt = result.x
# aircraft final design params


original_values = [
    original_aircraft['geo_param']['wing']['AR'],
    original_aircraft['geo_param']['wing']['S'],
    original_aircraft['geo_param']['wing']['sweep'],
    original_aircraft['geo_param']['EH']['Cht'],
    original_aircraft['dimensions']['ldg']['xnlg'],
    original_aircraft['dimensions']['ldg']['xmlg'],
    original_aircraft['dimensions']['ldg']['ymlg'],
    references['W0'],

]

optimized_values = [
    aircraft['geo_param']['wing']['AR'],
    aircraft['geo_param']['wing']['S'],
    aircraft['geo_param']['wing']['sweep'],
    aircraft['geo_param']['EH']['Cht'],
    aircraft['dimensions']['ldg']['xnlg'],
    aircraft['dimensions']['ldg']['xmlg'],
    aircraft['dimensions']['ldg']['ymlg'],
    result.fun * references['W0']
]

print(original_values)
print(optimized_values)

df = pd.DataFrame({
    'Design Variable': ['AR_w', 'Sw', 'sweep_w', 'Cht', 'xnlg', 'xmlg', 'ymlg', 'W0'],
    'Original Value': original_values,
    'Optimized Value': optimized_values,
    'Percent Change (%)': [100 * (opt - orig) / orig for orig, opt in zip(original_values, optimized_values)]
})

print(df)

df.to_excel('PRJ-23\\Otimização\\Resultados\\Potimizacao.xlsx', index=False)


# generating a 3d plot
fig_3d, ax_3d = plot3d([aircraft, original_aircraft], labels=['Optimized', 'Original'], show=True)
ax_3d.view_init(90, 0)
ax_3d.set_zticks([])
fig_3d.savefig('PRJ-23\\Otimização\\Resultados\\3dplot.png', dpi=500)


# Plot optimization history
plt.style.use('default')
plt.rcParams['font.family'] = 'Segoe UI'  # Set font to Segoe UI
#sns.set_palette("Set2")

# Plot 1: Design variables
fig1, ax1 = plt.subplots(figsize=(12, 6))
Xarray = np.array(Xlist)
design_var = [r'$AR_w$', r'$S_w$', r'$\Lambda_w$', r'$C_{ht}$', r'$x_{nlg}$', r'$x_{mlg}$', r'$y_{mlg}$']
for i in range(Xarray.shape[1]):
    ax1.plot(Xarray[:, i], label=design_var[i], linewidth=2)
ax1.set_ylabel('Variáveis de Design', fontsize=14)
ax1.set_xlabel('Iterações', fontsize=14)
ax1.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, shadow=True)
ax1.grid(True, alpha=0.3)
sns.despine()
plt.tight_layout()
plt.savefig('PRJ-23\\Otimização\\Resultados\\design_variables_history.png', dpi=500)

# Plot 2: Objective function
fig2, ax2 = plt.subplots(figsize=(12, 6))
ax2.plot(flist, color='red', linewidth=2)
ax2.set_ylabel(r'Função Objetivo ($W_0$)', fontsize=14)
ax2.set_xlabel('Iterações', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
sns.despine()
plt.tight_layout()
plt.savefig('PRJ-23\\Otimização\\Resultados\\objective_function_history.png', dpi=500)

# Plot 3: Constraints
fig3, ax3 = plt.subplots(figsize=(12, 8))
garray = np.array([g[0] for g in glist])  # Todos os itens exceto o último
constraint_labels = [
    r'$\Delta S_{wlan} \geq 0$',
    r'$SM_{fwd} \leq 0.30$',
    r'$SM_{aft} \geq 0.05$',
    r'$f_{nlg,fwd} \leq 0.18$',
    r'$f_{nlg,aft} \geq 0.03$',
    r'$\alpha_{tipback} \geq 15°$',
    r'$\alpha_{tailstrike} \geq 10°$',
    r'$\phi_{overturn} \leq 63°$',
    r'$b_{tank}/b_w \leq 0.95$',
    r'$b_w \leq 36m$',
    r'$y_{mlg} \leq 4.5m$',
    r'$z_{tv} - z_{mlg} \leq 13.5m$',
    r'$L_h \geq 13m$',
    # r'$|x_{EV} - x_{EH}| \leq 0.5m$'
]

for i in range(garray.shape[1]):
    ax3.plot(garray[:, i], label=constraint_labels[i], linewidth=2)
ax3.axhline(y=0, color='black', linestyle='--', alpha=0.7, linewidth=2)
ax3.set_ylabel('Restrições', fontsize=14)
ax3.set_xlabel('Iterações', fontsize=14)
ax3.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, shadow=True)
ax3.grid(True, alpha=0.3)
sns.despine()
plt.tight_layout()
plt.savefig('PRJ-23\\Otimização\\Resultados\\constraints_history.png', dpi=500)



print('Restrições finais desnormalizadas: ')
print(f'Delta S_wlan: {confun(Xopt)[0]* references["deltaS_wlan"]}')
pprint(confun(Xopt))



print("\n" + "="*30 + " TABELA DE VARIÁVEIS DE PROJETO " + "="*30)

# 1. Tabela de Variáveis (Batentes)
var_names_clean = ['AR_w', 'Sw', 'sweep_w', 'Cht', 'xnlg', 'xmlg', 'ymlg']
var_data = []
tol_vars = 1e-3

for i, val_norm in enumerate(Xopt):
    var_name = var_names_clean[i]
    
    # Recalcula o valor real
    if var_name == 'sweep_w':
        lower, upper = bounds[var_name]
        real_val = val_norm * (upper - lower) + lower
    else:
        lower, upper = bounds[var_name]
        real_val = val_norm * (upper - lower) + lower

    # Verifica Status
    if val_norm <= tol_vars:
        status = "BATENTE INF."
    elif val_norm >= (1.0 - tol_vars):
        status = "BATENTE SUP."
    else:
        status = "OK"
        
    var_data.append([var_name, status, round(val_norm, 4), round(real_val, 4), f"[{lower}, {upper}]"])

df_vars_final = pd.DataFrame(var_data, columns=['Variável', 'Status', 'Valor Norm.', 'Valor Real', 'Limites'])
print(df_vars_final.to_string(index=False))


print("\n" + "="*30 + " TABELA DE RESTRIÇÕES " + "="*30)

# 2. Tabela de Restrições (Ativas)
final_cons_values = confun(Xopt)
cons_data = []
tol_cons = 1e-4

for i, val in enumerate(final_cons_values):
    label = constraint_labels[i].replace('$', '').replace('\\', '') # Limpa LaTeX
    
    if val < -tol_cons:
        status = "VIOLADA"
    elif abs(val) <= tol_cons:
        status = "ATIVA"
    else:
        status = "INATIVA (Folga)"
        
    cons_data.append([label, status, round(val, 6)])

df_cons_final = pd.DataFrame(cons_data, columns=['Restrição', 'Status', 'Margem Final'])
print(df_cons_final.to_string(index=False))
print("\n")