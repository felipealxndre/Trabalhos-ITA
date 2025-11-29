#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pymoo.core.problem import Problem
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.sampling.lhs import LHS
from aux_tools_doe import corrdot
import design_tools_myaircraft as dt

#=========================================

# SETUP

def des_funcs(Xinp):

    # Fixed parameters
    # Ponto de Projeto
    g = 9.81 # Aceleração da gravidade
    T0_guess = 125600 #Chute inicial
    W0_guess = 490000.0 #Chute inicial
    altitude_cruise = 11000.0000 #Carteado
    Mach_cruise = 0.7700000 #Range de 0.75 a 0.80
    range_cruise = 3700e3 # Req projeto
    range_altcruise = 370400 # 200 NM
    loiter_time = 2700.00000 # 45 minutos
    altitude_altcruise = 4572.00000 # Caso de testess
    Mach_altcruise = 0.40000000 #Carteado
    distance_takeoff = 1800.0 #Req projeto
    distance_landing = 1150.0 # Req projeto

    #Parametros carteados
    altitude_takeoff = 0
    TO_flap_def = 20 * np.pi / 180
    TO_slat_def = 0
    altitude_landing = 0
    LD_flap_def = 40 * np.pi / 180
    LD_slat_def = 0
    MLW_frac = 0.84

    S_w, AR_w, sweep_w, delta_w, xr_w = Xinp
    aircraft = dt.my_aircraft()  # Define a aeronave padrão.
    aircraft['geo_param']['wing'].update({'S': S_w})
    aircraft['geo_param']['wing'].update({'AR': AR_w})
    aircraft['geo_param']['wing'].update({'sweep': sweep_w})
    aircraft['geo_param']['wing'].update({'delta': delta_w})
    aircraft['geo_param']['wing'].update({'xr': xr_w})

    W0, Wf,_,_,_,SM_aft,_,_,_,_,_,_ = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise,Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

    # Returns
    return W0, Wf, SM_aft

# Lower and upper bounds of each input variable
aircraft = dt.my_aircraft()  # Define a aeronave padrão.
S_w = aircraft['geo_param']['wing']['S']
xr_w = aircraft['geo_param']['wing']['xr']

#Wing area, AR, Sweep, Dihedral, Wing root position, respectively
lb = [S_w*0.9 , 6 , 10*np.pi/180, 0*np.pi/180, xr_w*0.9]
ub = [S_w*1.1 , 14, 40*np.pi/180, 5*np.pi/180, xr_w*1.1]

# Desired number of samples
n_samples = 40

# Sampling type
#sampler = FloatRandomSampling()
sampler = LHS()

# Plot type (0-simple, 1-complete)
plot_type = 1

#=========================================

# EXECUTION

# Set random seed to make results repeatable
np.random.seed(123)

# Get the number of input variables
n_var = len(lb)

# Initialize problem with lower and upper
problem = Problem(n_var=n_var, xl=lb, xu=ub)

# Draw samples
X = sampler(problem, n_samples).get("X")

# Execute all cases and store outputs
W0_samples = np.zeros(n_samples)
Wf_samples = np.zeros(n_samples)
SM_samples = np.zeros(n_samples)

for ii in range(n_samples):

    # Evaluate sample
    (W0,Wf,SM_aft) = des_funcs(X[ii,:])

    # Store the relevant information
    W0_samples[ii] = W0/1e3  # Convert to kN
    Wf_samples[ii] = Wf/1e3  # Convert to kN
    SM_samples[ii] = SM_aft

# Create a pandas dataframe with all the information
df = pd.DataFrame({'S_w' : X[:,0],
                   'AR_w' : X[:,1],
                   'sweep_w' : X[:,2]*180/np.pi, #Convert to degrees
                   'delta_w' : X[:,3]*180/np.pi, #Convert to degrees
                   'xr_w' : X[:,4],
                   'W0' : W0_samples,
                   'Wf' : Wf_samples,
                   'SM_aft' : SM_samples})

# Plot the correlation matrix
sns.set_style(style='white', font_scale=1.4)

if plot_type == 0:

    # Simple plot
    fig = sns.pairplot(df,corner=True)

elif plot_type == 1:

    # Complete plot
    fig = sns.PairGrid(df, diag_sharey=False)
    fig.map_lower(sns.regplot, lowess=True, line_kws={'color': 'black'})
    fig.map_diag(sns.histplot)
    fig.map_upper(corrdot)

# Plot window
plt.tight_layout()
plt.show()
 
fig.savefig('Aula Lab 1/Analyze/Correlation_chart_40s_myaircraft.svg')
