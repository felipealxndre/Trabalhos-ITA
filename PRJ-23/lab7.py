"""
Gerador de Diagrama V-n (Envelope de Voo)
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

# configuração padrão
aircraft = dt.default_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)


# colocar os dados da aeronave:
gravity = 9.80665 # mesma do auxtools  

rho = {
    'SL': atmosphere(z=0, Tba=288.15)[2],
    '10000': atmosphere(z=10000 * 0.3048, Tba=288.15)[2],
    '35000': atmosphere(z=35000 * 0.3048, Tba=288.15)[2],
}

S = aircraft['geo_param']['wing']['S']         # Área da asa (m²)
c = aircraft['dimensions']['wing']['cm']       # Corda média aerodinâmica 

CL_max = 1.30161                               # Coeficiente de sustentação máximo (CL,max) - AVL Lab 4

CL_min = -0.5                                  # Coeficiente de sustentação máximo negativo (CL,min)- Lab 4

a = 7.35940                                    # dCL/dalpha [em radianos] - Lab 4

# Esses aqui são dados que pegamos do FAR Part 25

n_limit_pos = 2.5                               # Fator de carga limite positivo 
                                                # O regulamento (FAR 25.337) define 2.5 para aeronaves de categoria transporte
n_limit_neg = -1.0                              # Fator de carga limite negativo (n_limit, neg). 
                                                # O regulamento (FAR 25.337) define -1.0 para aeronaves de transporte
VC = 221.36                                     # Velocidade de Cruzeiro de Projeto (VC) - antes tava 175.0 
                                                # A aeronave deve ser capaz de suportar rajadas nesta velocidade
VD = 1.25 * VC                                  # Velocidade de Mergulho de Projeto (VD). 
                                                # O regulamento exige VD ≥ 1.25 * VC

# Parâmetros de Rajada (Gust) - também vem do FAR 25

Ude_VC = 56 * 0.3048                            # Velocidade de Rajada de Projeto (Ude) em VC
                                                # 56 ft/s em VC do regulamento
Ude_VD = 28 * 0.3048                            # Velocidade de Rajada de Projeto (Ude) em VD
                                                # 28 ft/s em VD (metade de VC) do regulamento


def calcular_diagrama_vn(W, rho, nome_caso):
    """
    Calcula e plota o diagrama V-n para uma condição de voo
    
    W: Peso da aeronave em N
    rho: Densidade do ar em kg/m³
    nome_caso: Nome descritivo da condição
    """
    
    print(f"\n{'='*70}")
    print(f"  {nome_caso}")
    print(f"{'='*70}")
    
    W_S = W / S
    print(f"\n  Peso (W):              {W:,.0f} N")
    print(f"  Densidade (ρ):         {rho:.3f} kg/m³")
    print(f"  Carga Alar (W/S):      {W_S:.2f} N/m²")
    # Velocidade de Estol (1g) -> L = W
    VS = math.sqrt((2 * W_S) / (rho * CL_max))
    # Velocidade de Manobra
    VA = VS * math.sqrt(n_limit_pos)
    
    try:
        # Velocidade de Interseção do Estol Negativo com o Limite Estrutural
        V_neg_intersect = VS * math.sqrt(
            abs(n_limit_neg) / abs(CL_min / CL_max)
        )
        V_neg_intersect = min(V_neg_intersect, VD)
    except:
        V_neg_intersect = VD
    
    print(f"\n  VELOCIDADES:")
    print(f"    V_S (Estol, 1g):       {VS:.2f} m/s")
    print(f"    V_A (Manobra):         {VA:.2f} m/s")
    print(f"    V_C (Cruzeiro):        {VC:.2f} m/s")
    print(f"    V_D (Mergulho):        {VD:.2f} m/s")
    
    # valores de velocidade para a curva de estol positivo
    v_pos_stall = np.linspace(0, VA, 100)
    # n = (V / V_S) ** 2
    n_pos_stall = (v_pos_stall / VS) ** 2
    # valores de velocidade para a curva de estol negativo
    v_neg_stall = np.linspace(0, V_neg_intersect, 100)
    # n = (CL_min / CL_max) * (V / V_S) ** 2
    n_neg_stall = (CL_min / CL_max) * (v_neg_stall / VS) ** 2
    
    # limites estruturais
    maneuver_V = [VA, VD, VD, V_neg_intersect, 0]
    maneuver_n = [n_limit_pos, n_limit_pos, n_limit_neg, n_limit_neg, 
                  n_neg_stall[0] if len(n_neg_stall) > 0 else 0]
    # fator de massa e alívio
    mu_g = (2 * W_S) / (rho * c * a * gravity)
    Kg = (0.88 * mu_g) / (5.3 + mu_g)
    # fator de carga incremental em V_C e V_D
    delta_n_C = (Kg * Ude_VC * VC * a * rho) / (2 * W_S)
    delta_n_D = (Kg * Ude_VD * VD * a * rho) / (2 * W_S)
    
    print(f"\n  RAJADA:")
    print(f"    Fator de Massa (μ_g):  {mu_g:.3f}")
    print(f"    Fator de Alívio (K_g): {Kg:.3f}")
    print(f"    Δn em V_C:             {delta_n_C:+.3f}")
    print(f"    Δn em V_D:             {delta_n_D:+.3f}")
    print(f"{'='*70}\n")
    # linhas de rajada
    gust_V = [0, VC, VD]
    gust_n_pos = [1, 1 + delta_n_C, 1 + delta_n_D]
    gust_n_neg = [1, 1 - delta_n_C, 1 - delta_n_D]
    
    plt.figure(figsize=(14, 9))
    
    plt.plot(v_pos_stall, n_pos_stall, 'b-', linewidth=2, label='Estol Positivo ($C_{L,max}$)')
    plt.plot(v_neg_stall, n_neg_stall, 'b--', linewidth=2, label='Estol Negativo ($C_{L,min}$)')
    plt.plot(maneuver_V, maneuver_n, 'r-', linewidth=2.5, label='Limite Estrutural')
    
    plt.plot(gust_V, gust_n_pos, 'g-.', linewidth=2, 
            label=f'Rajada Positiva ({Ude_VC:.1f} m/s @ Vc)')
    plt.plot(gust_V, gust_n_neg, 'g-.', linewidth=2, 
            label=f'Rajada Negativa ({Ude_VC:.1f} m/s @ Vc)')
    
    plt.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    plt.axhline(y=0, color='black', linewidth=0.5)
    
    plt.axvline(x=VA, color='purple', linestyle=':', alpha=0.5)
    plt.text(VA, n_limit_neg - 0.2, f'$V_A$ = {VA:.1f} m/s', 
            ha='center', fontsize=10, color='purple')
    
    plt.text(VC, n_limit_neg - 0.2, f'$V_C$ = {VC:.1f} m/s', 
            ha='center', fontsize=10, color='green')
    
    plt.text(VD, n_limit_neg - 0.2, f'$V_D$ = {VD:.1f} m/s', 
            ha='center', fontsize=10, color='red')
    
    plt.title(f'Diagrama V-n Combinado\n{nome_caso}', fontsize=16, fontweight='bold')
    plt.xlabel('Velocidade Equivalente, $V_e$ (m/s)', fontsize=13)
    plt.ylabel('Fator de Carga, $n$', fontsize=13)
    plt.xlim(0, VD + 25)
    plt.ylim(n_limit_neg - 0.5, n_limit_pos + 0.5)
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'Resultados/diagrama_vn_{nome_caso.replace(" ", "_").lower()}.png', dpi=300)
    
    # Pegar o dataframe do caso
    df_caso = get_dataframe(nome_caso, VS, VA, VC, VD, W, rho, mu_g, Kg, delta_n_C, delta_n_D)

    return df_caso

def get_dataframe(nome_caso, VS, VA, VC, VD, W, rho, mu_g, Kg, delta_n_C, delta_n_D):
    
    # Criar DataFrame com uma linha por caso
    df = pd.DataFrame({
        'Caso': [nome_caso],
        'Peso_N': [W],
        'Densidade_kg_m3': [rho],
        'Carga_Alar_N_m2': [W/S],
        'VS_ms': [VS],
        'VA_ms': [VA],
        'VC_ms': [VC],
        'VD_ms': [VD],
        'Mu_g': [mu_g],
        'Kg': [Kg],
        'Delta_n_C': [delta_n_C],
        'Delta_n_D': [delta_n_D]
    })
    
    return df


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
# W0, Wf, We, deltaS_wlan, SM_fwd, SM_aft, b_tank_b_w, frac_nlg_fwd, frac_nlg_aft, alpha_tipback, alpha_tailstrike, phi_overturn = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)
pesos = dt.analyze(aircraft, W0_guess, T0_guess, Mach_cruise, altitude_cruise, range_cruise, Mach_altcruise, range_altcruise, altitude_altcruise, loiter_time, altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def, altitude_landing, distance_landing, LD_flap_def, LD_slat_def, MLW_frac)

MTOW = pesos[0]  # peso máximo de decolagem
MZFW = pesos[0] - pesos[1]  # peso máximo sem combustível = MTOW - peso do combustível
W_projeto = 39071.40*gravity  # ponto de projeto - cruzeiro típico

print(f"MTOW: {MTOW:.2f} N")
print(f"MZFW: {MZFW:.2f} N")
print(f"W_projeto: {W_projeto:.2f} N")

# casos de estudo
# dai temos que saber os pesos e densidades em cada caso
caso1 = calcular_diagrama_vn(
    W=MTOW,
    rho=rho['SL'],
    nome_caso="MTOW @ Nível do Mar"
)

caso2 = calcular_diagrama_vn(
    W=MZFW,
    rho=rho['SL'],
    nome_caso="MZFW @ Nível do Mar"
)

caso3 = calcular_diagrama_vn(
    W=W_projeto,
    rho=rho['35000'],
    nome_caso="Peso de Cruzeiro @ 35000 ft"
)

caso4 = calcular_diagrama_vn(
    W=W_projeto,
    rho=rho['10000'],
    nome_caso="Peso de Pouso @ 10000 ft"
)

df = pd.concat([caso1, caso2, caso3, caso4], ignore_index=True)
df.to_excel('Resultados/Diagrama_Vn_Completo.xlsx', index=False)