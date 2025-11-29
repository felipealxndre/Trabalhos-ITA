import numpy as np

# --- DADOS EXTRAÍDOS DOS CARTÕES DE VOO ---
# Voo 1 (CG Dianteiro)
t_voo1_yaw_dumper_on = np.array([0, 3.1])
n_picos_voo1_yaw_dumper_on = 3

t_voo1_yaw_dumper_off = np.array([0, 13])
n_picos_voo1_yaw_dumper_off = 10

# Voo 2 (CG Traseiro)
t_voo2_yaw_dumper_on = np.array([0, 2.81])
n_picos_voo2_yaw_dumper_on = 3

t_voo2_yaw_dumper_off = np.array([0, 9.09])
n_picos_voo2_yaw_dumper_off = 7

# --- CÁLCULOS ---


def analisar_dinamica(tempos, n_picos):

    n_ciclos = max(1, n_picos - 1)

    Td = 2*(tempos[1] - tempos[0])/n_ciclos

    wd = (2 * np.pi) / Td

    n_picos = np.ceil(n_picos/2)
    zeta = (7 - n_picos) / 10.0

    wn = wd / np.sqrt(1 - zeta**2)

    return Td, wn, zeta


Td1_yaw_dumper_on, wn1_yaw_dumper_on, zeta1_yaw_dumper_on = analisar_dinamica(
    t_voo1_yaw_dumper_on, n_picos_voo1_yaw_dumper_on)
Td1_yaw_dumper_off, wn1_yaw_dumper_off, zeta1_yaw_dumper_off = analisar_dinamica(
    t_voo1_yaw_dumper_off, n_picos_voo1_yaw_dumper_off)
Td2_yaw_dumper_on, wn2_yaw_dumper_on, zeta2_yaw_dumper_on = analisar_dinamica(
    t_voo2_yaw_dumper_on, n_picos_voo2_yaw_dumper_on)
Td2_yaw_dumper_off, wn2_yaw_dumper_off, zeta2_yaw_dumper_off = analisar_dinamica(
    t_voo2_yaw_dumper_off, n_picos_voo2_yaw_dumper_off)

print("-" * 40)
print("ANÁLISE EXPERIMENTAL - DUTCH ROLL")
print("-" * 40)
print(f"VOO 1 (CG Dianteiro):")
print("-" * 40)
print("Yaw Dumper ON")
print(f"  Tempos entre picos: {Td1_yaw_dumper_on:.2f} s")
print(f"  Frequência Natural (wn): {wn1_yaw_dumper_on:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1_yaw_dumper_on:.4f}")
print("-" * 40)
print("Yaw Dumper OFF")
print(f"  Tempos entre picos: {Td1_yaw_dumper_off:.2f} s")
print(f"  Frequência Natural (wn): {wn1_yaw_dumper_off:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1_yaw_dumper_off:.4f}")
print("-" * 40)
print(f"VOO 2 (CG Traseiro):")
print("-" * 40)
print("Yaw Dumper ON")
print(f"  Tempos entre picos: {Td2_yaw_dumper_on:.2f} s")
print(f"  Frequência Natural (wn): {wn2_yaw_dumper_on:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2_yaw_dumper_on:.4f}")
print("-" * 40)
print("Yaw Dumper OFF")
print(f"  Tempos entre picos: {Td2_yaw_dumper_off:.2f} s")
print(f"  Frequência Natural (wn): {wn2_yaw_dumper_off:.4f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2_yaw_dumper_off:.4f}")
print("-" * 40)
