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
import aux_tools as at
import seaborn as sns

# configuração padrão
aircraft = dt.default_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)


# colocar os dados da aeronave:
gravity = 9.80665   
RHO_SL = 1.225     # Densidade nível do mar (kg/m³)

S = 120.77         # Área da asa (m²)
c = 5.64           # Corda média aerodinâmica (m)
CL_max = 2.1       # CLmax positivo
CL_min = -0.8      # CLmax negativo
a = 5.9            # É o CLalpha (1/rad)

n_limit_pos = 2.5  # Fator de carga limite positivo
n_limit_neg = -1.0 # Fator de carga limite negativo

VC = 175.0         # Velocidade de Cruzeiro (m/s)
VD = 218.5         # Velocidade de Mergulho (m/s)

Ude_VC = 17.07     # Velocidade da rajada de projeto em V_C (m/s)
Ude_VD = 8.53      # Velocidade da rajada de projeto em V_D (m/s)


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
    # Velocidade de Estol (1g) -- L = W
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

# casos de estudo
# dai temos que saber os pesos e densidades em cada caso
caso1 = calcular_diagrama_vn(
    W=800544.0,
    rho=RHO_SL,
    nome_caso="MTOW @ Nível do Mar"
)


caso2 = calcular_diagrama_vn(
    W=667120.0,
    rho=RHO_SL,
    nome_caso="MZFW @ Nível do Mar"
)

caso3 = calcular_diagrama_vn(
    W=733616.0,
    rho=0.381,
    nome_caso="Peso de Cruzeiro @ 35000 ft"
)

caso4 = calcular_diagrama_vn(
    W=711712.0,
    rho=1.056,
    nome_caso="Peso de Pouso @ 10000 ft"
)

df = pd.concat([caso1, caso2, caso3, caso4], ignore_index=True)
df.to_excel('Resultados/Diagrama_Vn_Completo.xlsx', index=False)