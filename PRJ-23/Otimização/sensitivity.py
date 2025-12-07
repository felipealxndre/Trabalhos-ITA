#Grupo - Pinguins de Madagascar
#Integrante:
# - Thiago Akira Missato

# Importar os arquivos
import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt
from scipy.stats import linregress

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
import design_tools_optmized as dt
import aux_tools as at

aircraft = dt.my_aircraft() #Define a aeronave nossa.
new_dimensions = dt.geometry(aircraft) #Calcula as dimensões da aeronave.
aircraft['dimensions'].update(new_dimensions)   #Atualiza as dimensões da aeronave.

#Caso de testes



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

res_ref = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise,
                     Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, 
                     altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, 
                     altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

# Desempacotando referência e criando dicionário
ref_keys = ['W0', 'Wf', 'We', 'deltaS_wlan', 'SM_fwd', 'SM_aft', 'b_tank_b_w', 
            'frac_nlg_fwd', 'frac_nlg_aft', 'alpha_tipback', 'alpha_tailstrike', 'phi_overturn']
references = dict(zip(ref_keys, res_ref))

# =============================================================================
# 2. DEFINIÇÃO DOS PARÂMETROS DE ESTUDO
# =============================================================================

# Captura os valores iniciais para a lista de iteração
# vals = {
#     'Sw': aircraft['geo_param']['wing']['S'],
#     'AR_w': aircraft['geo_param']['wing']['AR'],
#     'sweep_w': aircraft['geo_param']['wing']['sweep'],
#     'delta_w': aircraft['geo_param']['wing']['dihedral'],
#     'xr_w': aircraft['geo_param']['wing']['xr'],
#     'tcr_w': aircraft['geo_param']['wing']['tcr'],
#     'tct_w': aircraft['geo_param']['wing']['tct'],
#     'Cht': aircraft['geo_param']['EH']['Cht'],
#     'Lc': aircraft['geo_param']['EH']['Lc'],
#     'Cvt': aircraft['geo_param']['EV']['Cvt'],
#     'Lb': aircraft['geo_param']['EV']['Lb'],
#     'x_tailstrike': aircraft['data']['misc']['x_tailstrike'],
#     'z_tailstrike': aircraft['data']['misc']['z_tailstrike'],
#     'xnlg': aircraft['dimensions']['ldg']['xnlg'],
#     'xmlg': aircraft['dimensions']['ldg']['xmlg'],
#     'ymlg': aircraft['dimensions']['ldg']['ymlg'],
#     'Mach_cruise': Mach_cruise,
#     'range_cruise': range_cruise
# }

vals = {
    'Sw': aircraft['geo_param']['wing']['S'],
    'AR_w': aircraft['geo_param']['wing']['AR'],
    'sweep_w': aircraft['geo_param']['wing']['sweep'],
    'Cht': aircraft['geo_param']['EH']['Cht'],
    'xnlg': aircraft['dimensions']['ldg']['xnlg'],
    'xmlg': aircraft['dimensions']['ldg']['xmlg'],
    'ymlg': aircraft['dimensions']['ldg']['ymlg']
}

param_names = list(vals.keys())
input_vals = list(vals.values())
output = []

print(f"\n--- Iniciando Análise de Sensibilidade para {len(param_names)} variáveis ---\n")


for i, name in enumerate(param_names):
    # Reseta a aeronave para o estado limpo a cada iteração
    aircraft = dt.my_aircraft()
    
    # Valores padrão de operação para resetar
    current_mach = Mach_cruise
    current_range = range_cruise
    
    # Perturbação padrão de 2%
    perturbation = 1.02 
    perc_change = 0.02 # Usado para dividir (delta Y / delta X)
    
    # Aplicação da perturbação
    if name == 'Sw':
        aircraft['geo_param']['wing']['S'] *= perturbation
    elif name == 'AR_w':
        aircraft['geo_param']['wing']['AR'] *= perturbation
    elif name == 'sweep_w':
        aircraft['geo_param']['wing']['sweep'] *= perturbation
    elif name == 'delta_w':
        # Para diedro, se for 0, adicionamos um valor absoluto pequeno (ex: 2 graus)
        # Se não for 0, multiplicamos
        val_atual = aircraft['geo_param']['wing']['dihedral']
        if val_atual == 0:
            aircraft['geo_param']['wing']['dihedral'] = 2 * np.pi / 180
            # Evita divisão por zero: aqui a sensibilidade não é %/%, é absoluta
            perc_change = 1.0 # Placeholder, requer análise manual se diedro for 0
        else:
            aircraft['geo_param']['wing']['dihedral'] *= perturbation
    elif name == 'xr_w':
        aircraft['geo_param']['wing']['xr'] *= perturbation
    elif name == 'tcr_w':
        aircraft['geo_param']['wing']['tcr'] *= perturbation
    elif name == 'tct_w':
        aircraft['geo_param']['wing']['tct'] *= perturbation
    elif name == 'Cht':
        aircraft['geo_param']['EH']['Cht'] *= perturbation
    elif name == 'Lc':
        aircraft['geo_param']['EH']['Lc'] *= perturbation
    elif name == 'Cvt':
        aircraft['geo_param']['EV']['Cvt'] *= perturbation
    elif name == 'Lb':
        aircraft['geo_param']['EV']['Lb'] *= perturbation
    elif name == 'x_tailstrike':
        aircraft['data']['misc']['x_tailstrike'] *= perturbation
    elif name == 'z_tailstrike':
        aircraft['data']['misc']['z_tailstrike'] *= perturbation
    elif name == 'xnlg':
        aircraft['dimensions']['ldg']['xnlg'] *= perturbation
    elif name == 'xmlg':
        aircraft['dimensions']['ldg']['xmlg'] *= perturbation
    elif name == 'ymlg':
        aircraft['dimensions']['ldg']['ymlg'] *= perturbation
    elif name == 'Mach_cruise':
        current_mach *= perturbation
    elif name == 'range_cruise':
        current_range *= perturbation

    # Recalcula geometria com a perturbação
    new_dimensions = dt.geometry(aircraft)
    aircraft['dimensions'].update(new_dimensions)
    
    # Roda análise
    res_new = dt.analyze(aircraft, W0_guess, T0_guess, current_mach, altitude_cruise, current_range,
                         Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, 
                         altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, 
                         altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
    
    # Calcula sensibilidade Relativa: (% variação Saída) / (% variação Entrada)
    current_output = []
    
    for k, key in enumerate(ref_keys):
        val_ref = references[key]
        val_new = res_new[k]
        
        # Proteção contra divisão por zero na referência (ex: margem estática pode ser 0)
        if abs(val_ref) < 1e-6:
            # Se a referência é zero, calculamos variação absoluta simples
            sens = (val_new - val_ref)
        else:
            # Sensibilidade adimensional: (dY/Y) / (dX/X)
            pct_diff_y = (val_new - val_ref) / val_ref
            sens = pct_diff_y / perc_change
            
        current_output.append(round(sens, 4))

    print(f"Calculado: {name}")
    output.append(current_output)


df_sensitivity = pd.DataFrame(output, index=param_names, columns=ref_keys)

# Exibe no console

print(df_sensitivity[['W0', 'Wf', 'deltaS_wlan', 'SM_fwd', 'SM_aft']].to_string()) # Mostra só as principais pra conferir

# Salva arquivo
df_sensitivity.to_excel('sensibilidade.xlsx')