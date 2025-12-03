import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt
from scipy.stats import linregress

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import design_tools as dt
import aux_tools as at

aircraft = dt.my_aircraft() #Define a aeronave padrão.

aircraft['geo_param']['wing']['AR'] = 8.585856
aircraft['geo_param']['wing']['S'] = 82.913206
aircraft['geo_param']['wing']['sweep'] = 0.349066
aircraft['geo_param']['EH']['Cht'] = 1.300000 
aircraft['dimensions']['ldg']['xnlg'] = 2.808354 
aircraft['dimensions']['ldg']['xmlg'] = 16.832371
aircraft['dimensions']['ldg']['ymlg'] = 2.470000 


new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)

#Ponto de Projeto

gravity = 9.81
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

TO_flap_def = 20 * np.pi / 180 
LD_flap_def = aircraft['data']['flap']['max_def']
TO_slat_def = 0
LD_slat_def = 0
h_ground = 10.668
altitude_cruise = 11000
altitude_takeoff = 0.0
altitude_landing = 0.0
MLW_frac = 0.84

# Pista de decolagem e pouso
distance_takeoff = 1800.0 #Req projeto
distance_landing = 1150.0 # Req projeto

#Parametros carteados
TO_flap_def = 20 * np.pi / 180 
LD_flap_def = aircraft['data']['flap']['max_def']
TO_slat_def = 0
LD_slat_def = 0
h_ground = 35.0*0.3048
altitude_cruise = 11000
altitude_takeoff = 0.0
altitude_landing = 0.0
MLW_frac = 0.84

W0, We, Wf, xcg_e, T0, T0vec, S_wlan = dt.thrust_matching(aircraft, W0_guess, T0_guess,
                                                           TO_flap_def, LD_flap_def,
                                                           TO_slat_def, LD_slat_def,
                                                           h_ground,
                                                           altitude_cruise, Mach_cruise, range_cruise,
                                                           loiter_time,
                                                           altitude_altcruise, Mach_altcruise, range_altcruise,
                                                           altitude_takeoff, distance_takeoff,
                                                           altitude_landing, distance_landing, MLW_frac)

S_req = S_wlan # Área da asa de pouso necessária
deltaS_wlan = aircraft['geo_param']['wing']['S'] - S_wlan
print(f'Delta S_wlan = {deltaS_wlan} m²') # Diferença entre a área da asa atual e a necessária
print(f'S_wlan = {S_req} m²') # Área da asa de pouso necessária

# Variação da área da asa 

Sw_range = np.linspace(90.0, 130.0, 100) # Faixa de valores para Sw
Swlan_range = np.zeros(len(Sw_range)) # Inicializa o vetor de área da asa de pouso
#wing_charge = np.zeros(len(Sw_range)) # Inicializa o vetor de carga alar
T0vec_range = np.zeros((len(Sw_range),len(T0vec))) # Inicializa o vetor de T0
i = 0

for Sw in Sw_range:
    aircraft['geo_param']['wing'].update({'S': Sw}) # Atualiza a área da asa
    new_dimensions = dt.geometry(aircraft)  #Calcula as dimensões da aeronave.
    aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave. 
    W0, We, Wf, xcg_e, T0, T0vec, S_wlan = dt.thrust_matching(aircraft, W0_guess, T0_guess,
                                                           TO_flap_def, LD_flap_def,
                                                           TO_slat_def, LD_slat_def,
                                                           h_ground,
                                                           altitude_cruise, Mach_cruise, range_cruise,
                                                           loiter_time,
                                                           altitude_altcruise, Mach_altcruise, range_altcruise,
                                                           altitude_takeoff, distance_takeoff,
                                                           altitude_landing, distance_landing, MLW_frac)
    T0vec_range[i,:] = T0vec # Armazena o valor de T0
    #wing_charge[i] = W0/Sw  # Calcula a carga alar
    Swlan_range[i] = S_wlan  # Armazena a área da asa de pouso
    i = i+1



labels = ['Takeoff','Cruise','FAR25.111','FAR25.121a','FAR25.121b','FAR25.121c','FAR25.119','FAR25.121d']
plt.figure(figsize=(16,10))  # Aumenta o tamanho da figura
T0vec_range = T0vec_range/1e3  # Converter para kN
plt.plot(Sw_range, T0vec_range, '-', linewidth=3)

vert_y = np.linspace(58,136, 100)
vert_x = np.zeros(len(vert_y))
vert_x = vert_x + 1*S_req  # Define a linha vertical em S_req
plt.plot(vert_x, vert_y, '-', linewidth=3, label='Área de asa de pouso necessária') 

# Adiciona ponto vermelho na coordenada (100, tração de cruzeiro correspondente)
area_ponto = S_req*1.05
print(f'Área de asa de pouso necessária com margem: {area_ponto} m²')  # Área de asa de pouso necessária com margem
# Encontra o índice mais próximo de 100 em Sw_range
idx_ponto = np.argmin(np.abs(Sw_range - area_ponto))
tracao_cruise = T0vec_range[idx_ponto, 0]  # índice 1 = 'Cruise' conforme labels
print(f'Tração de cruzeiro com margem: {tracao_cruise*1.05} kN')  # Tração de cruzeiro em kN
plt.plot(area_ponto, tracao_cruise*1.05, 'ro', markersize=10, label='Ponto de projeto')

# Calcula o maior valor entre as duas curvas em cada ponto
upper = np.maximum(T0vec_range[:, 0], T0vec_range[:, 4])

# Limita a hachura também pela linha vertical em S_req
mask = (Sw_range >= S_req)

plt.fill_between(
    Sw_range[mask],
    upper[mask],
    T0vec_range.max() * 1.05,
    color='none',
    hatch='xxx',
    edgecolor='red',
    linewidth=0.1,
    alpha=0.4,
    zorder=0
)

plt.xlabel('Área de asa (m²)', fontsize=22)
plt.ylabel('Tração Total(kN)', fontsize=22)
plt.tick_params(axis='both', labelsize=22)
plt.legend(
    labels + ['Landing'] + ['Ponto de Projeto'],
    loc='lower center',
    bbox_to_anchor=(0.5, 1.02),  # centraliza e coloca fora, acima do gráfico
    fontsize=18,
    ncol=5,
)

plt.grid(True)
plt.savefig('PRJ-23\\Resultados\\prjdiagram.png')

