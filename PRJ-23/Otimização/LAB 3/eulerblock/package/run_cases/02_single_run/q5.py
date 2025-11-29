'''
This script is an example of how to use the EulerBlock interface
to analyze an airfoil.

Written by:
Ney Rafael Secco
Instituto Tecnologico de Aeronautica
Sao Jose dos Campos - Brazil
ney@ita.br
Sept 2022
'''

# Modify path to include the Euler solver
import sys
with open('../eulerblock_path.txt') as f:
    exec(f.read())
sys.path.append(eulerblock_path)
from eulerblock import euler_mod as eb
import airfoil_mod as am

import numpy as np
import matplotlib.pyplot as plt

# Flow conditions
alpha_1 = 0.9403848319587892*np.pi/180.0
alpha_2 = 1.186816407 * np.pi/180.0
alphas = [alpha_1, alpha_2]
alpha_labels = ['Aerofólio Otimizado', 'NACA1411']
colors = ['blue', 'red']
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

# Airfoil CST values - ClarkY
#Al = [-0.1294, -0.0036, -0.0666]
#Au = [0.206, 0.2728, 0.2292]

# Airfoil CST values - NACA1411
Al_naca = [-0.1489439, -0.10330027, -0.10305128, -0.10514982]
Au_naca = [ 0.16146332,  0.18349204,  0.14126241,  0.18194397]

# Optimized airfoil CST variables
Al_opt = [-0.05, -0.20314275, -0.10647595, -0.05]
Au_opt = [0.12694442, 0.13648835, 0.18591317, 0.18153676]

# Airfoil CST values - NACA0012
#Al = [-0.16941, -0.15145, -0.13904, -0.13988]
#Au = [ 0.16941,  0.15145,  0.13904,  0.13988]

# Airfoil CST values - RAE2822 (remember to lower alpha)
# Au = [0.12569744, 0.14564925, 0.1520378,  0.21346977, 0.17867926, 0.20872853]
# Al = [-0.13334359, -0.11827498, -0.22179445, -0.12644608, -0.08233977,  0.05193736]

# Store results for both cases
results_list = []
airfoils_list = []

# Run for both angles of attack
for i, alpha in enumerate(alphas):
    print(f'\nRunning case {i+1}: alpha = {alpha*180/np.pi:.3f} degrees')

    if alpha == alpha_2:
        Al = Al_naca
        Au = Au_naca
    else:
        Al = Al_opt
        Au = Au_opt
    
    results = eb.run_cst(Al, Au, Nchord, alpha, mach,
                         gamma, order,
                         iter, dt, CFL, use_local_dt, res_NK, res_tol,
                         reinitialize, adj_funcs, plot=False, NJ=NJ, s0=s0)
    
    results_list.append(results)
    
    # Refine coordinates for each case
    x = (1-np.cos(np.linspace(0, 1, 81)*np.pi))/2
    airfoil = am.cstfoil(Au, Al, x)
    airfoils_list.append(airfoil)
    
    # Print results for this case
    print(f'Results for {alpha_labels[i]}:')
    print('max thickness=', results['maxt'])
    print('max thickness x/c=', results['xmaxt'])
    print('max camber=', results['maxc'])
    print('max camber x/c=', results['xmaxc'])
    print('CL=', results['CL'])
    print('CD=', results['CD'])
    print('CM=', results['CM'])
    print('grads=', results['grads'])

# Export airfoil (same for both cases)
am.export_airfoil(airfoils_list[0], title='AIRFOIL', filename='airfoil.dat')

# Plot overlayed curves
fig = plt.figure(figsize=(10, 8))

# Cp distribution
plt.subplot(311)
for i, (results, label, color) in enumerate(zip(results_list, alpha_labels, colors)):
    xx = results['distrib']['xx']
    Cp = results['distrib']['Cp']
    plt.plot(xx, Cp, color=color, label=label, linewidth=2)
plt.ylabel('Cp')
plt.legend()
plt.grid(True, alpha=0.3)
plt.gca().invert_yaxis()

# Mach distribution
plt.subplot(312)
for i, (results, label, color) in enumerate(zip(results_list, alpha_labels, colors)):
    xx = results['distrib']['xx']
    MachD = results['distrib']['Mach']
    plt.plot(xx, MachD, color=color, label=label, linewidth=2)
plt.ylabel('Mach')
plt.legend()
plt.grid(True, alpha=0.3)

# Airfoil coordinates (same for both cases)
plt.subplot(313)
for i, (airfoil, label, color) in enumerate(zip(airfoils_list, alpha_labels, colors)):
    xf = airfoil['x_coord']
    yf = airfoil['y_coord']
    plt.plot(xf, yf, color=color, label=label, linewidth=2)
plt.ylabel('y/c')
plt.xlabel('x/c')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# Summary table
print('\n' + '='*60)
print('SUMMARY COMPARISON')
print('='*60)
print(f"{'Parameter':<15} {'α₁ = 0.94°':<15} {'α₂ = 1.19°':<15}")
print('-'*45)
for param in ['CL', 'CD', 'CM']:
    val1 = results_list[0][param]
    val2 = results_list[1][param]
    print(f"{param:<15} {val1:<15.6f} {val2:<15.6f}")