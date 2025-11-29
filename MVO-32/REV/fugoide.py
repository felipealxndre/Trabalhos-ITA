import numpy as np

# --- DADOS EXTRAÍDOS DOS CARTÕES DE VOO ---

# Voo 1 (CG Dianteiro ~20.15%) - Fonte: Pág 11 do Cartão
t_voo1 = np.array([0, 25.1, 50.9, 60+11.1, 60+43.5, 120+13.3, 60+37.3])
Vi_voo1 = np.array([165, 200, 170, 195, 175, 190, 180])
Zpi_voo1 = np.array([17980, 17200, 17780, 17180, 17620, 17220, 17480])

# Voo 2 (CG Traseiro ~22.99%) - Fonte: Pág 13 do Cartão
t_voo2 = np.array([0, 21.9, 45.1, 60+10.05, 60+37.3, 120+6.5, 120+28.3])
Vi_voo2 = np.array([190, 170, 185, 175, 182, 176, 180])
Zpi_voo2 = np.array([15250, 15500, 15200, 15450, 15300, 15400, 15300])


def analisar_dinamica(tempos, amplitudes, n_picos):
    # Cálculo do decremento logarítmico
    A1 = abs(amplitudes[1] - amplitudes.mean())  # Primeira amplitude
    An = abs(amplitudes[-2] - amplitudes.mean())  # Última amplitude
    delta = (1/n_picos) * np.log(A1/An)

    # Cálculo do coeficiente de amortecimento
    zeta = 1/np.sqrt(1 + (2*np.pi/delta)**2)

    # Cálculo do período experimental
    T = (tempos[-2] - tempos[1])/n_picos

    # Cálculo da frequência natural não amortecida
    wn = 2*np.pi/(T*np.sqrt(1-zeta**2))

    return T, wn, zeta


T1, wn1, zeta1 = analisar_dinamica(t_voo1, Vi_voo1, n_picos=2)
T2, wn2, zeta2 = analisar_dinamica(t_voo2, Vi_voo2, n_picos=2)

print("-" * 40)
print("ANÁLISE EXPERIMENTAL - FUGÓIDE")
print("-" * 40)
print(f"VOO 1 (CG Dianteiro):")
print("-" * 40)
print(f"  Tempos entre picos: {T1:.2f} s")
print(f"  Frequência Natural (wn): {wn1:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1:.4f}")
print("-" * 40)
print(f"VOO 2 (CG Traseiro):")
print("-" * 40)
print(f"  Tempos entre picos: {T2:.2f} s")
print(f"  Frequência Natural (wn): {wn2:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2:.4f}")
print("-" * 40)
