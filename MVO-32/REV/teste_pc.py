import numpy as np

# --- DADOS EXTRAÍDOS DOS CARTÕES DE VOO ---

# Voo 1 (CG Dianteiro ~20.15%) - Fonte: Pág 11 do Cartão
t_voo1_pulsos = np.array([0, 2.2])
n_picos_voo1_pulsos = 3

t_voo1_doublet1 = np.array([0, 1])
n_picos_voo1_doublet1 = 2

t_voo1_doublet2 = np.array([0, 0.65])
n_picos_voo1_doublet2 = 2

# Voo 2 (CG Traseiro ~22.99%) - Fonte: Pág 13 do Cartão
t_voo2_pulsos = np.array([0, 1.94])
n_picos_voo2_pulsos = 3

t_voo2_doublet = np.array([0, 3.58])
n_picos_voo2_doublet = 3


def analisar_dinamica(tempos, n_picos):

    n_ciclos = max(1, n_picos - 1)

    Td = 2*(tempos[1] - tempos[0])/n_ciclos

    wd = (2 * np.pi) / Td

    n_picos = np.floor(n_picos/2)
    zeta = (7 - n_picos) / 10.0

    wn = wd / np.sqrt(1 - zeta**2)

    return Td, wn, zeta


Td1_pulsos, wn1_pulsos, zeta1_pulsos = analisar_dinamica(
    t_voo1_pulsos, n_picos_voo1_pulsos)
Td1_doublet1, wn1_doublet1, zeta1_doublet1 = analisar_dinamica(
    t_voo1_doublet1, n_picos_voo1_doublet1)
Td1_doublet2, wn1_doublet2, zeta1_doublet2 = analisar_dinamica(
    t_voo1_doublet2, n_picos_voo1_doublet2)
Td2_pulsos, wn2_pulsos, zeta2_pulsos = analisar_dinamica(
    t_voo2_pulsos, n_picos_voo2_pulsos)
Td2_doublet, wn2_doublet, zeta2_doublet = analisar_dinamica(
    t_voo2_doublet, n_picos_voo2_doublet)

print("-" * 40)
print("ANÁLISE EXPERIMENTAL - CURTO PERÍODO")
print("-" * 40)
print(f"VOO 1 (CG Dianteiro):")
print("-" * 40)
print("Pulsos")
print(f"  Tempos entre picos: {Td1_pulsos:.2f} s")
print(f"  Frequência Natural (wn): {wn1_pulsos:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1_pulsos:.4f}")
print("-" * 40)
print("Doublet 1")
print(f"  Tempos entre picos: {Td1_doublet1:.2f} s")
print(f"  Frequência Natural (wn): {wn1_doublet1:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1_doublet1:.4f}")
print("-" * 40)
print("Doublet 2")
print(f"  Tempos entre picos: {Td1_doublet2:.2f} s")
print(f"  Frequência Natural (wn): {wn1_doublet2:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1_doublet2:.4f}")
print("-" * 40)
print(f"VOO 2 (CG Traseiro):")
print("-" * 40)
print("Pulsos")
print(f"  Tempos entre picos: {Td2_pulsos:.2f} s")
print(f"  Frequência Natural (wn): {wn2_pulsos:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2_pulsos:.4f}")
print("-" * 40)
print("Doublet")
print(f"  Tempos entre picos: {Td2_doublet:.2f} s")
print(f"  Frequência Natural (wn): {wn2_doublet:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2_doublet:.4f}")
print("-" * 40)
