'''
INSTITUTO TECNOLÓGICO DE AERONÁUTICA
PROGRAMA DE ESPECIALIZAÇÃO EM ENGENHARIA AERONÁUTICA
OTIMIZAÇÃO MULTIDISCIPLINAR

Código com aplicações de ferramentas do Python para otimização
multi-objetivo (pymoo 0.6.0)

Maj. Ney Sêcco 2025
'''

# IMPORTS
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pprint import pprint

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import auxmod as am
import design_tools as dt
import aux_tools as at
from plot3d import plot3d

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

# EXECUTION

# Define optimization problem
class MyProblem(ElementwiseProblem):

    def __init__(self):

        # Set general characteristics of the problem
        super().__init__(n_var=7, # Number of design variables - 7 variáveis
                         n_obj=2, # Number of objective functions - W0, Wf
                         n_constr=13, # Number of constraints
                         # bounds normalizados [0,1] para todas as 7 variáveis
                         xl=[0, 0, 0, 0, 0, 0, 0], # Lower bounds of design variables
                         xu=[1, 1, 1, 1, 1, 1, 1]) # Upper bound of design variables

    def _evaluate(self, X, out, *args, **kwargs):

        # design variables - Sw, ARw, sweep, Cht, x_nlg, x_mlg, y_mlg
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

        # Compute objectives - W0 e Wf
        new_dimensions = dt.geometry(aircraft) #Calcula as dimensões da aeronave.
        aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave
        # Compute objective function

        # analyse to see objective functions
        W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
        f1 = W0
        f2 = Wf

        # Compute constraints - aqui tem que ser menor igual a zero
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
        ]

        # truque pra colocar todas como negativas pra ser no formato do pymoo
        constraints = [-c for c in constraints]

        # Gather results
        out["F"] = [f1, f2]
        out["G"] = constraints

# Create an instance of the Problem
problem = MyProblem()

# Select optimization algorithm - MOGA
# tamanho da população: pop size
algorithm = NSGA2(pop_size=200, eliminate_duplicates=True)


# Define reference values
aircraft = dt.my_aircraft() # Define a aeronave nossa como referencia.
new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave.

original_aircraft = aircraft.copy()

# creating the references - from the reference group aircraft
reference_keys = ['W0', 'Wf', 'We', 'deltaS_wlan', 'SM_fwd', 'SM_aft', 'b_tank_b_w', 'frac_nlg_fwd', 'frac_nlg_aft', 'alpha_tipback', 'alpha_tailstrike', 'phi_overturn']

reference_values = dt.analyze(
    aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise,
    Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time,
    altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def,
    altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac
)
references = dict(zip(reference_keys, reference_values))

# Define design of the 7 variable bounds - dictonary of lists - [lower, upper]
bounds = {
    'AR_w': [7, 12],
    'Sw': [80, 120],
    'sweep_w': [0, 40], # degrees
    'Cht': [1.3, 1.6],
    'xnlg': [1, 3],
    'xmlg': [15, 17],
    'ymlg': [2, 6]
}

# Solve the optimization
res = minimize(problem,
               algorithm,
               ('n_gen', 200), # número de gerações - quantas vezes ele vai replicar
               seed=1,
               verbose=True)

plt.style.use('default')
plt.rcParams['font.family'] = 'Segoe UI'  # Set font to Segoe UI
sns.set_palette("Set2")

plt.figure(figsize=(12, 5))
sns.scatterplot(x=res.F[:,0], y=res.F[:,1], s=50)
sns.lineplot(x=res.F[:,0], y=res.F[:,1])

# Sort the Pareto front by W0 (first objective) to get proper beginning, middle, end
sorted_indices = np.argsort(res.F[:, 0])
n_points = len(sorted_indices)
print(n_points)

# Select beginning (min W0), middle, and end (max W0) points
indices = [sorted_indices[0], sorted_indices[n_points//2], sorted_indices[-1]]
pareto_points = ['Ponto A', 'Ponto B', 'Ponto C']

colors = ['purple','red', 'orange']
for i, idx in enumerate(indices):
    plt.scatter(res.F[idx, 0], res.F[idx, 1], color=colors[i], s=50, zorder=5, marker= 'D', label=pareto_points[i])

plt.xlabel(r'$W_0$', fontsize=14)
plt.ylabel(r'$W_f$', fontsize=14)
sns.despine()
plt.tight_layout()
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig('PRJ-23\\Otimização\\Resultados\\pareto_front.png', dpi=500)

# Sort the Pareto front by W0 (first objective) to get proper beginning, middle, end
sorted_indices = np.argsort(res.F[:, 0])
n_points = len(sorted_indices)

# Select beginning (min W0), middle, and end (max W0) points
indices = [sorted_indices[0], sorted_indices[n_points//2], sorted_indices[-1]]
pareto_points = ['Ponto A', 'Ponto B', 'Ponto C']

# Create aircraft configurations for the 3 selected points
aircraft_configs = []
pareto_data = []

for i, idx in enumerate(indices):
    # Get design variables for this point
    X_point = res.X[idx]
    
    # Denormalize design variables
    AR_w = X_point[0] * (bounds['AR_w'][1] - bounds['AR_w'][0]) + bounds['AR_w'][0]
    S_w = X_point[1] * (bounds['Sw'][1] - bounds['Sw'][0]) + bounds['Sw'][0]
    sweep_w = (X_point[2] * (bounds['sweep_w'][1] - bounds['sweep_w'][0]) + bounds['sweep_w'][0]) * np.pi/180
    Cht = X_point[3] * (bounds['Cht'][1] - bounds['Cht'][0]) + bounds['Cht'][0]
    xnlg = X_point[4] * (bounds['xnlg'][1] - bounds['xnlg'][0]) + bounds['xnlg'][0]
    xmlg = X_point[5] * (bounds['xmlg'][1] - bounds['xmlg'][0]) + bounds['xmlg'][0]
    ymlg = X_point[6] * (bounds['ymlg'][1] - bounds['ymlg'][0]) + bounds['ymlg'][0]
    
    # Create aircraft copy and update parameters
    aircraft_pareto = aircraft.copy()
    aircraft_pareto['geo_param']['wing']['AR'] = AR_w
    aircraft_pareto['geo_param']['wing']['S'] = S_w
    aircraft_pareto['geo_param']['wing']['sweep'] = sweep_w
    aircraft_pareto['geo_param']['EH']['Cht'] = Cht
    aircraft_pareto['dimensions']['ldg']['xnlg'] = xnlg
    aircraft_pareto['dimensions']['ldg']['xmlg'] = xmlg
    aircraft_pareto['dimensions']['ldg']['ymlg'] = ymlg
    
    # Update geometry
    new_dimensions = dt.geometry(aircraft_pareto)
    aircraft_pareto['dimensions'].update(new_dimensions)
    
    aircraft_configs.append(aircraft_pareto)
    
    # Add data to DataFrame
    pareto_data.append({
        'Point': pareto_points[i],
        'AR_w': AR_w,
        'S_w': S_w,
        'sweep_w': sweep_w,
        'C_ht': Cht,
        'x_nlg': xnlg,
        'x_mlg': xmlg,
        'y_mlg': ymlg,
        'W0': res.F[idx, 0],
        'Wf': res.F[idx, 1]
    })
    
    print(f"{pareto_points[i]} Point - W0: {res.F[idx, 0]:.0f}, Wf: {res.F[idx, 1]:.0f}")


pareto = pd.DataFrame(pareto_data)
print("\nPareto Points Comparison:")
print(pareto.round(3))

# Save to Excel
pareto.to_excel('PRJ-23\\Otimização\\Resultadospareto_points_comparison.xlsx', index=False)

# Plot 3D views for the 3 selected points side by side, each compared to original
fig = plt.figure(figsize=(18, 6))

for i, (aircraft_config, point_name) in enumerate(zip(aircraft_configs, pareto_points)):
    ax = fig.add_subplot(1, 3, i+1, projection='3d')
    
    # Plot both original and current configuration
    fig_temp, ax_temp = plot3d([aircraft_config, original_aircraft], 
                               labels=[point_name, 'Original'], show=False)

    # Copy the plot data to the subplot
    for line in ax_temp.lines:
        xdata = line._verts3d[0]
        ydata = line._verts3d[1] 
        zdata = line._verts3d[2]
        ax.plot(xdata, ydata, zdata, color=line.get_color(), label=line.get_label())
    
    ax.set_title(f'{point_name} vs Original', fontsize=14)
    ax.view_init(90, 0)
    ax.set_zticks([])
    ax.legend()
    
    plt.close(fig_temp)  # Close temporary figure

plt.tight_layout()
plt.savefig('PRJ-23\\Otimização\\Resultados\\pareto_3d_comparison.png', dpi=500)