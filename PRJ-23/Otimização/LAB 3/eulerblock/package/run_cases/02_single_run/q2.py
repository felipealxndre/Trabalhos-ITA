import sys
with open('../eulerblock_path.txt') as f:
    exec(f.read())
sys.path.append(eulerblock_path)
from eulerblock import euler_mod as eb
import airfoil_mod as am

import numpy as np
import matplotlib.pyplot as plt

# Flow conditions
alphas = np.linspace(-5, 5, 5)*np.pi/180.0  # Angle of attack [rad]
mach = 0.75 # mach number
gamma = 1.4 # ratio of specific heats of the fluid

# Solver parameters
order = 2 # 1 for 1st; 2 for 2nd
iter = 20000 # Use 50000 for first order
dt = 0.001 # Global time step (only used if use_local_dt==0) .Use 0.002 for first order and 0.001 for second order.
CFL = 0.2 # Courant number used to compute local time steps (only used if use_local_dt==1)
use_local_dt = 1 # Time stepping-strategy. 0 for global dt; 1 for local dt based on CFL.
res_NK = 1e-4
res_tol = 1e-6
reinitialize = 0
adj_funcs = ['cl_jlow', 'cd_jlow', 'cm_jlow']

# Discretization (higher grid levels have higher number of elements)
grid_level = 1.0

# Reference mesh properties
Nchord0 = 30 # Number of elements along the chord
NJ0 = 48 # Number of elements along the normal direction
s00 = 0.5e-2 # Reference mesh spacing

Nchord = int((Nchord0*grid_level) + 1)
NJ = int((NJ0*grid_level) + 1)
s0 = NJ0/(NJ-1)*s00

# Airfoil CST values - NACA1411
Al = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]
Au = [ 0.16146332,  0.18349204,  0.14126241,  0.18194397]

CLs = []

for alpha in alphas:
    results = eb.run_cst(Al, Au, Nchord, alpha, mach,
                         gamma, order,
                         iter, dt, CFL, use_local_dt, res_NK, res_tol,
                         reinitialize, adj_funcs, plot=False, NJ=NJ, s0=s0)

    CLs.append(results['CL'])


import numpy as np
import seaborn as sns
plt.figure()
sns.scatterplot(x = alphas*180.0/np.pi, y = CLs)
# linear interpolating
slope, intercept = np.polyfit(alphas*180.0/np.pi, CLs, 1)
# plotting the line
sns.lineplot(x=alphas*180.0/np.pi, y=intercept + slope * (alphas*180.0/np.pi), color='red', label = r'$C_L = $' + f'{slope:.4f}' + r'$\cdot\alpha + $' + f'{intercept:.4f}')

plt.xlabel('Ângulo de Ataque (graus)')
plt.ylabel('Coeficiente de Sustentação CL')
plt.show()