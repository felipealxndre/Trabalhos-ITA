# IMPORTS
import numpy as np
from scipy.optimize import minimize
import auxmod as am
import matplotlib.pyplot as plt
import design_tools_default as dt
from pprint import pprint
from plot3d import plot3d
import copy  # Add this import
# removing warnings in terminal
import warnings
import pandas as pd
warnings.filterwarnings("ignore")

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

# Define history list
Xlist = []
flist = []
g1list = []
g2list = []


# Define objective function
def objfun(X):

    # design variables - Sw, ARw, sweep, Cht, X_nlg, X_mlg, Y_mlg

    # desnormalizando
    AR_w = X[0] * (bounds['AR_w'][1] - bounds['AR_w'][0]) + bounds['AR_w'][0]
    S_w = X[1] * (bounds['Sw'][1] - bounds['Sw'][0]) + bounds['Sw'][0]

    # inputing other restrictions
    aircraft['geo_param']['wing']['AR'] = AR_w
    aircraft['geo_param']['wing']['S'] = S_w

    new_dimensions = dt.geometry(aircraft) #Calcula as dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave
    # Compute objective function

    # analyse to see objective functions
    W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

    # Não acho necessário normalizar a função objetivo;
    f = W0

    # Save results for the plot
    Xlist.append(X)
    flist.append(f)

    return f

# Define constraints
def confun(X):

    # design variables - Sw, ARw, sweep, Cht, xnlg, xmlg, ymlg

    # desnormalizando
    AR_w = X[0] * (bounds['AR_w'][1] - bounds['AR_w'][0]) + bounds['AR_w'][0]
    S_w = X[1] * (bounds['Sw'][1] - bounds['Sw'][0]) + bounds['Sw'][0]

    # inputing other restrictions
    aircraft['geo_param']['wing']['AR'] = AR_w
    aircraft['geo_param']['wing']['S'] = S_w

    new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave
    
    # Compute objective function

    # analyse to see objective functions
    W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
    
    # Constraints
    constraints = [
        (30 - aircraft['dimensions']['wing']['b']) / 30,        # b_w <= 30 m
    ]
    # Optionally, save results for plotting
    g1list.append(constraints[0])

    return constraints

# Create list of constraints

con1 = {
    'type': 'ineq',
    'fun': confun
}

cons = [con1]

# Define starting point and reference values 
aircraft = dt.my_aircraft() # Define a aeronave nossa como referencia.
new_dimensions = dt.geometry(aircraft) # Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   # Atualiza as dimensões da aeronave.

original_aircraft = copy.deepcopy(aircraft)  # Use deepcopy instead of copy

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
}

# normalizados -> value = (value - lower) /(upper - lower)
X0 = {
    'AR_w': (7.5 - bounds['AR_w'][0])/(bounds['AR_w'][1] - bounds['AR_w'][0]),
    'Sw': (90 - bounds['Sw'][0])/(bounds['Sw'][1] - bounds['Sw'][0]),
}

X0 = np.array(list(X0.values()))
print('Starting Point: \n')
print(X0)



# Additional options
options = {'maxiter':200, 'ftol':1e-6}

# Run optimizer
result = minimize(objfun, X0,
                  constraints=cons, 
                  bounds= [[0, 1]] * 2, # normalized bounds - all between 0 ad 1
                  method='slsqp', options=options)

# Print results
print(result)
Xopt = result.x
# aircraft final design params


original_values = [
    7.5,
    90,
    np.sqrt(90 * 7.5),
]

optimized_values = [
    aircraft['geo_param']['wing']['AR'],
    aircraft['geo_param']['wing']['S'],
    aircraft['dimensions']['wing']['b'],
]

print(original_values)
print(optimized_values)

df = pd.DataFrame({
    'Design Variable': ['AR_w', 'Sw', 'b_w'],
    'Original Value': original_values,
    'Optimized Value': optimized_values,
    'Percent Change (%)': [100 * (opt - orig) / orig for orig, opt in zip(original_values, optimized_values)]
})

print(df)

df.to_excel('otimizacao_parte1.xlsx', index=False)

# generating a 3d plot
plot3d([aircraft, original_aircraft], labels=['Optimized', 'Original'], show=True)
