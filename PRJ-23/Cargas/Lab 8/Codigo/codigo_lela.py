"""
Laboratório de Cargas em Aeronaves
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import design_tools as dt
from aux_tools import atmosphere
import seaborn as sns
import sympy as sp


# configuração padrão
aircraft = dt.my_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)


# colocar os dados da aeronave:
gravity = 9.80665 # mesma do auxtools  

rho = {
    'SL': atmosphere(z=0, Tba=288.15)[2],
    '10000': atmosphere(z=10000 * 0.3048, Tba=288.15)[2],
    '35000': atmosphere(z=35000 * 0.3048, Tba=288.15)[2],
}

T = {
    'SL': atmosphere(z=0, Tba=288.15)[0],
    '10000': atmosphere(z=10000 * 0.3048, Tba=288.15)[0],
    '35000': atmosphere(z=35000 * 0.3048, Tba=288.15)[0],
}

#Dados de referência
# calculando os pesos

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

n_limit_pos = 2.5

# W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
#pesos = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

MTOW = 254000  # peso máximo de decolagem (kg)
MZFW = 192700  # peso máximo sem combustível = MTOW - peso do combustível

VMO = 340.0  # Velocidade máxima operacional (m/s) KEAS
MMO = 0.9    # Mach máximo operacional
v_s1g = 165.0  # Velocidade de estol com 1g (m/s) KEAS

W = MTOW*gravity
P_limit = n_limit_pos*W #carga limite - talvez tenha que trocar MTOW para W para poder testar vários casos
P_ultimate = P_limit*1.5 #carga final

print(f"W = {W:.2f} N")
print(f"P_limit = {P_limit:.2f} N")
print(f"P_ultimate = {P_ultimate:.2f} N")

VA = 134.2 #m/s
rho_sl = rho['SL']
#para a nossa aeronave:
## características aerodinâmicas e geométricas da aeronave
# S = aircraft['geo_param']['wing']['S']         # Área da asa (m²)
# c = aircraft['dimensions']['wing']['cm']       # Corda média aerodinâmica 
# Cm_acw = -0.15  #vem do AVL!!!- CORRIGIR # Momento em torno do CA (acw - aerodynamic center wing)

# #posições importantes
# xcg = aircraft['weights']['xcg_list'][5]
# xac_w = aircraft['weights']['xac_w']
# xac_h = aircraft['weights']['xac_h']
## características aerodinâmicas e geométricas da aeronave
S = 360        # Área da asa (m²)
c = 8.3        # Corda média aerodinâmica 
Cm_acw = -0.15  #vem do AVL!!!- CORRIGIR # Momento em torno do CA (acw - aerodynamic center wing)

#posições importantes
xcg = 38
xac_w = 39
xac_h = 68
lw = xac_w-xcg
lh = xac_h - xcg
#equilíbrio de forças

q = 0.5 * rho_sl * VA**2
Mac_w = q*S*c*Cm_acw

print(f"q = {q:.2f} Pa")
print(f"Mac_w = {Mac_w:.2f} N.m")
print(f"lw = {lw:.2f} m")
print(f"lh = {lh:.2f} m")

A = np.array([[1.0, 1.0],
              [lw,  lh]])
b = np.array([n_limit_pos*W, -Mac_w])  # note o sinal: passa Mac_w para o outro lado
Lw, Lh = np.linalg.solve(A, b)

print(f"Lw = {Lw:.5f} MN")
print(f"Lh = {Lh:.5f} MN")

#parte 2
MTOW = 79000  # peso máximo de decolagem (kg)
W = MTOW*gravity
W_wfsemi = 0.23*W  # peso na semi-asa
Weng = 2400 *gravity  # peso do motor
yeng = 5.10  # posição do motor
n_limit_pos = 2.5
VA = 245*0.514444 #m/s
Ltotal = n_limit_pos*W
Lsemi = Ltotal/2
print(f"W = {W:.2f} N")
print(f"Ltotal = {Ltotal:.2f} N")
print(f"Lsemi = {Lsemi:.2f} N")

#geometria da aeronave
b2 = 17.15   # semi-envergadura [m]
S = 125       # área da asa [m²]
cr = 6.07    # corda raiz [m]
lam = 0.2    # afilamento
ct = lam*cr  # corda ponta

# variável ao longo da semi-asa
y = sp.symbols('y', real=True)

# ----- Distribuição de corda (linear) -----
c = cr + (ct - cr)*y/b2  # decresce linearmente
# constantes que escalam a distribuição (proporcional à c(y))

k_L = Lsemi/sp.integrate(c, (y, 0, b2))   # V(y)
k_w = W_wfsemi/sp.integrate(c, (y, 0, b2))   # V(y)

L_prime = k_L * c    # carga aerodinâmica
W_prime = k_w * c
# carga líquida
N_prime = L_prime - W_prime
eta = sp.symbols('eta', real=True)
N_prime_eta = N_prime.subs(y, eta)
#y>=yeng
V1 = sp.integrate(N_prime_eta, (eta, y, b2))   # V(y)
sp.simplify(V1)
#y<yeng
V2 = sp.integrate(N_prime_eta, (eta, y, b2)) + Weng   # V(y)
sp.simplify(V2)

print(f"Distribuição de carga V1 (y >= {yeng} m):")
sp.pprint(sp.simplify(V1))
print(f"\nDistribuição de carga V2 (y < {yeng} m):")
sp.pprint(sp.simplify(V2))
print(f"\nValores chaves para V(y):")
print(f"V({b2}) = {V1.subs(y, b2):.2f} N/m")
print(f"V({yeng}+) = {V1.subs(y, yeng):.2f} N/m")
print(f"V({yeng}-) = {V2.subs(y, yeng):.2f} N/m") #conferir
print(f"V(0) = {V2.subs(y, 0):.2f} N/m")
    
#calculo de momento fletor
zeta = sp.symbols('zeta', real=True)
V1_zeta = V1.subs(y, zeta)
V2_zeta = V2.subs(y, zeta)
M1 = sp.integrate(V1_zeta, (zeta, y, b2)) #y>=yeng
M2 = sp.integrate(V1_zeta, (zeta, yeng, b2)) + sp.integrate(V2_zeta, (zeta, y, yeng)) #y<yeng

print(f"\nValores chaves para M(y):") 
print(f"M({b2}) = {M1.subs(y, b2):.2f} N.m")

#Momento distribuído
N_prime_y = N_prime*y
Mdist = sp.integrate(N_prime_y, (y, 0, b2))
print(f"\nMomento distribuído na semi-asa: M_dist = {Mdist:.2f} N.m")

#Momento do motor
Meng = Weng * yeng

#Momento total na raiz 
Mtotal = Mdist + Meng
print(f"\nMomento total na raiz da semi-asa: M_total = {Mtotal:.2f} N.m")   

sp.plotting.plot((M1, (y, 0, yeng)),(M2,(y, yeng, b2)), title='Momento Fletor na Semi-Asa', xlabel='y (m)', ylabel='M(y) (N.m)', legend=False)   
sp.plotting.plot((V2, (y, 0, yeng)), (V1, (y, yeng, b2)), title='Esforço cortante na Semi-Asa', xlabel='y (m)', ylabel='V(y) (N.m)', legend=False)
