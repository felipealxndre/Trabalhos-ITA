import math
import numpy as np
import matplotlib.pyplot as plt
import design_tools as dt
from aux_tools import atmosphere

# conversão de unidades
MPS_TO_KT = 1.94384  # 1 m/s ≈ 1.94384 nós

# configuração padrão
aircraft = dt.my_aircraft()
new_dimensions = dt.geometry(aircraft)
aircraft['dimensions'].update(new_dimensions)

n_limit_pos = 2.5
n_limit_neg = -1.0

g = 9.80665 

# (rho, T) para cada condição
condicoes_voo = {
    'Sea Level': (atmosphere(z=0, Tba=288.15)[2], 288.15),
    '10000':     (atmosphere(z=10000 * 0.3048, Tba=268.34)[2], 268.34),
    '35000':     (atmosphere(z=35000 * 0.3048, Tba=218.8)[2], 218.8),
}

S = aircraft['geo_param']['wing']['S']         # Área da asa (m²)
c = aircraft['dimensions']['wing']['cm']       # Corda média aerodinâmica 

CL_max = 1.30161                               # Coeficiente de sustentação máximo (CL,max) - AVL Lab 4
CL_min = -0.5                                  # Coeficiente de sustentação máximo negativo (CL,min)- Lab 4
a      = 7.35940                               # dCL/dalpha [em radianos] - Lab 4


def calcular_Ude(W):
    """
    Modelo Embraer (Equivalent Height):

      heq = (W/S) / 10700
      Ude_VC = 56 ft/s * exp(-heq)
      Ude_VD = 0.5 * Ude_VC

    Retorna Ude_VC e Ude_VD em m/s.
    """
    W_S = W / S  # N/m²

    heq = W_S / 10700.0

    Ude_VC_fps = 56.0 * math.exp(-heq)   # ft/s
    Ude_VD_fps = 0.5  * Ude_VC_fps       # ft/s

    U_de_VC = Ude_VC_fps * 0.3048        # m/s
    U_de_VD = Ude_VD_fps * 0.3048        # m/s

    return U_de_VC, U_de_VD


# Velocidades de Projeto (em TAS)
def velocidades_projeto(W, rho_ar, temp_ar):
    W_S = W / S
    
    # Velocidade de Estol (TAS)
    V_s = math.sqrt((2 * W_S) / (rho_ar * CL_max))
    n_s = 0.5 * rho_ar * V_s**2 * CL_max / W_S

    # Velocidade de Manobra (TAS)
    V_A = V_s * math.sqrt(n_limit_pos)
    n_A = n_limit_pos
    
    # Velocidade de Cruzeiro (TAS)
    V_C = Mach_cruise * math.sqrt(1.4 * 287.0 * temp_ar)

    # Velocidade de Mergulho (TAS equivalente à VD em EAS)
    V_D_min = 1.25 * V_C 
    
    return V_s, n_s, V_A, n_A, V_C, V_D_min


def pico_rajada_1_cosseno(W_S, rho_ar, V_e, U_de):
    """
    Calcula o pico de fator de carga n para um dado V_e e U_de
    usando perfil de rajada 1-cosseno:

      U(s) = 0.5 * Ude * (1 - cos(pi * s / H)),  H = 25 * c

    Inclui fator de alívio Kg.
    """
    # Fator de alívio de rajada
    mu_g = (2 * W_S) / (rho_ar * c * a * g)
    K_g  = (0.88 * mu_g) / (5.3 + mu_g)

    # Densidade ao nível do mar para EAS
    rho_SL = condicoes_voo["Sea Level"][0]
    term_common = (rho_SL * K_g * a) / (2 * W_S)

    # Perfil 1-cosseno no espaço
    H = 25.0 * c
    s = np.linspace(0.0, 2.0 * H, 400)  # percorre a rajada inteira
    U_s = 0.5 * U_de * (1.0 - np.cos(np.pi * s / H))

    # n(s) = 1 + Δn(s) com Δn(s) = term_common * V_e * U(s)
    n_s = 1.0 + term_common * V_e * U_s

    return float(np.max(n_s))


def envelope_rajada(W, rho_ar, V_eC, V_eD):
    """
    Calcula n_C e n_D (rajada positiva) em função de VE,
    usando Ude com Equivalent Height + perfil 1-cosseno.
    """
    W_S = W / S
    
    # Ude em m/s via Equivalent Height
    U_de_VC, U_de_VD = calcular_Ude(W)

    # Pico de n em VC e VD com 1-cosseno
    n_C = pico_rajada_1_cosseno(W_S, rho_ar, V_eC, U_de_VC)
    n_D = pico_rajada_1_cosseno(W_S, rho_ar, V_eD, U_de_VD)
    
    return n_C, n_D


def plot_envelope_combinado(titulo,
                            V_es, V_eA, V_eC, V_eD,
                            n_C, n_D):
    """
    Plota o diagrama V-n combinado + envelope em preto
    seguindo a sequência geométrica desejada.
    """

    plt.figure(figsize=(14, 9))

    # 1) Vs- (velocidade onde estol negativo atinge n_limit_neg, limitada por VD)
    try:
        V_neg_intersect = V_es * math.sqrt(abs(n_limit_neg) / abs(CL_min / CL_max))
        V_neg_intersect = min(V_neg_intersect, V_eD)
    except:
        V_neg_intersect = V_eD
    
    # 2) Curvas de estol em VE (para visualização)
    v_pos_stall = np.linspace(0, V_eA, 200)
    n_pos_stall = (v_pos_stall / V_es) ** 2
    
    v_neg_stall = np.linspace(0, V_neg_intersect, 200)
    n_neg_stall = (CL_min / CL_max) * (v_neg_stall / V_es) ** 2
    
    plt.plot(v_pos_stall, n_pos_stall, 'b-.', linewidth=2, label='Estol')
    plt.plot(v_neg_stall, n_neg_stall, 'b-.', linewidth=2)

    # 3) Limites estruturais (manobra)
    plt.hlines(n_limit_pos, V_eA, V_eD, colors='red', linestyles='-.',
               linewidth=2.5, label='Limite Estrutural')
    plt.hlines(n_limit_neg, V_neg_intersect, V_eD, colors='red',
               linestyles='-.', linewidth=2.5)
    plt.vlines(V_eD, n_limit_neg, n_limit_pos, colors='red',
               linestyles='-.', linewidth=2.5)
    
    # 4) Envelope de rajada (linhas verdes principais – ligando 0–VC–VD)
    V_gust = np.array([0.0, V_eC, V_eD])
    n_gust_pos = np.array([1.0, n_C, n_D])
    n_gust_neg = np.array([1.0, 2.0 - n_C, 2.0 - n_D])

    plt.plot(V_gust, n_gust_pos, "g-.", linewidth=2, label='Rajada')
    plt.plot(V_gust, n_gust_neg, "g-.", linewidth=2)

    # Linhas tracejadas verdes retas de (0,1) até (VD, n_D) e (VD, 2-n_D)
    V_gust_full = np.linspace(0, V_eD, 200)
    n_gust_pos_line = 1 + (n_D - 1) * (V_gust_full / V_eD)
    n_gust_neg_line = 1 + ((2.0 - n_D) - 1) * (V_gust_full / V_eD)

    plt.plot(V_gust_full, n_gust_pos_line, "g-.", linewidth=2, alpha=0.8)
    plt.plot(V_gust_full, n_gust_neg_line, "g-.", linewidth=2, alpha=0.8)

    # Linhas de referência
    plt.axhline(1.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.axhline(0, color="black", linewidth=0.5)

    # Marcações de Vs, VA, VC e VD (linhas verticais + texto)
    plt.axvline(V_es, color="purple", linestyle=":", alpha=0.5)
    plt.axvline(V_eA, color="purple", linestyle=":", alpha=0.5)
    plt.axvline(V_eC, color="purple", linestyle=":", alpha=0.5)
    plt.axvline(V_eD, color="purple", linestyle=":", alpha=0.5)
    
    plt.text(V_es, n_limit_neg, f'$V_S$ = {V_es:.1f} m/s',
             ha="center", fontsize=12, color="purple",
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
    plt.text(V_eA, n_limit_neg, f'$V_A$ = {V_eA:.1f} m/s',
             ha="center", fontsize=12, color="purple",
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
    plt.text(V_eC, n_limit_neg, f'$V_C$ = {V_eC:.1f} m/s',
             ha="center", fontsize=12, color="purple",
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
    plt.text(V_eD, n_limit_neg, f'$V_D$ = {V_eD:.1f} m/s',
             ha="center", fontsize=12, color="purple",
             bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))

    # -----------------------------
    # ENVELOPE EM PRETO
    # -----------------------------

    # Grid de velocidade para o envelope
    V_env = np.linspace(V_es, V_eD, 400)

    # Rajadas interpoladas (usadas para recorte do envelope)
    gust_pos_interp = np.interp(V_env, V_gust, n_gust_pos)
    gust_neg_interp = np.interp(V_env, V_gust, n_gust_neg)

    n_top = np.zeros_like(V_env)
    n_bottom = np.zeros_like(V_env)

    for i, V in enumerate(V_env):
        # Estol positivo
        n_stall_pos = (V / V_es)**2

        # TOPO DO ENVELOPE
        if V <= V_eA:
            n_top[i] = n_stall_pos
        elif V <= V_eC:
            n_top[i] = max(n_limit_pos, gust_pos_interp[i])
        else:
            n_top[i] = max(n_limit_pos, gust_pos_interp[i])

        # FUNDO DO ENVELOPE
        if V <= V_neg_intersect:
            n_bottom[i] = 0.0
        elif V <= V_eC:
            n_bottom[i] = min(n_limit_neg, gust_neg_interp[i])
        else:
            n_struct = n_limit_neg + (0.0 - n_limit_neg) * (V - V_eC) / (V_eD - V_eC)
            n_bottom[i] = min(n_struct, gust_neg_interp[i])

    # curvas contínuas do envelope
    plt.plot(V_env, n_top,    color="black", linewidth=2,
             label='Envelope V-n Combinado')
    plt.plot(V_env, n_bottom, color="black", linewidth=2)

    # curvas verticais em Vs, Vs- e VD para fechar o contorno visualmente

    # vertical em Vs: de 0 até 1
    plt.plot([V_es, V_es], [0.0, 1.0], color="black", linewidth=2)

    # vertical em Vs-: de 0 até n_limit_neg
    plt.plot([V_neg_intersect, V_neg_intersect],
             [0.0, n_limit_neg], color="black", linewidth=2)

    # vertical em VD: de n_bottom(VD) até n_top(VD)
    plt.plot([V_eD, V_eD],
             [n_bottom[-1], n_top[-1]], color="black", linewidth=2)

    plt.xlabel("Velocidade Equivalente, $V_e$ (m/s)", fontsize=13)
    plt.ylabel("Fator de Carga, $n$", fontsize=13)
    plt.title(f'Diagrama V-n Combinado\n{titulo}', fontsize=16, fontweight='bold')
    plt.xlim(0, V_eD * 1.1)
    # plt.ylim(n_limit_neg - 0.5, n_limit_pos + 0.5)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best', fontsize=10)

    plt.tight_layout()
    plt.savefig(f"Diagrama_vn_{titulo.replace(' ', '_').lower()}.png", dpi=300)


def rodar_caso(W, rho_key, titulo):
    """
    Calcula velocidades (TAS), converte para VE (EAS) conforme a densidade,
    imprime os pontos principais em m/s e kt e plota o envelope.
    """
    rho_ar, temp_ar = condicoes_voo[rho_key]

    # velocidades de projeto em TAS
    V_s, n_s, V_A, n_A, V_C_local, VD_min = velocidades_projeto(W, rho_ar, temp_ar)

    # fator de conversão TAS -> VE
    rho_SL = condicoes_voo['Sea Level'][0]
    fator = math.sqrt(rho_ar / rho_SL)

    # em SL, fator = 1; em altitude, converte para VE
    V_es = V_s * fator
    V_eA = V_A * fator
    V_eC = V_C_local * fator
    V_eD = VD_min * fator

    # rajada em função de VE (via Ude + 1-cosseno)
    n_C, n_D = envelope_rajada(W, rho_ar, V_eC, V_eD)

    # conversão para nós
    V_es_kt = V_es * MPS_TO_KT
    V_eA_kt = V_eA * MPS_TO_KT
    V_eC_kt = V_eC * MPS_TO_KT
    V_eD_kt = V_eD * MPS_TO_KT

    # print dos pontos principais
    print(f"\n=== {titulo} ===")
    print(f"  rho  = {rho_ar:.4f} kg/m³")
    print(f"  V_S  = {V_es:.2f} m/s  ({V_es_kt:.1f} kt)")
    print(f"  V_A  = {V_eA:.2f} m/s  ({V_eA_kt:.1f} kt)")
    print(f"  V_C  = {V_eC:.2f} m/s  ({V_eC_kt:.1f} kt)")
    print(f"  V_D  = {V_eD:.2f} m/s  ({V_eD_kt:.1f} kt)")
    print(f"  n_C  = {n_C:.3f} (rajada em V_C)")
    print(f"  n_D  = {n_D:.3f} (rajada em V_D)")
    print(f"  n_C- = {2.0 - n_C:.3f} (rajada negativa em V_C)")
    print(f"  n_D- = {2.0 - n_D:.3f} (rajada negativa em V_D)")

    # plota para esse caso
    plot_envelope_combinado(
        titulo,
        V_es, V_eA, V_eC, V_eD,
        n_C, n_D
    )


# ----------------------------------------
# Cálculo de pesos
# ----------------------------------------

W0_guess = 43090 * g
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

pesos = dt.analyze(
    aircraft, W0_guess, T0_guess,
    Mach_cruise, altitude_cruise, range_cruise,
    Mach_altcruise, range_altcruise, altitude_altcruise,
    loiter_time,
    altitude_takeoff, distance_takeoff, TO_flap_def, TO_slat_def,
    altitude_landing, distance_landing, LD_flap_def, LD_slat_def,
    MLW_frac
)

MTOW      = pesos[0]            # peso máximo de decolagem
MZFW      = pesos[0] - pesos[1] # peso máximo sem combustível
W_projeto = 39071.40 * g       # ponto de projeto - cruzeiro típico


# ----------------------------------------
# Casos
# ----------------------------------------

rodar_caso(MTOW,      'Sea Level', "MTOW @ Nível do Mar")
rodar_caso(MZFW,      'Sea Level', "MZFW @ Nível do Mar")
rodar_caso(W_projeto, '10000',     "Peso de Cruzeiro @ 10000 ft")
rodar_caso(W_projeto, '35000',     "Peso de Cruzeiro @ 35000 ft")
