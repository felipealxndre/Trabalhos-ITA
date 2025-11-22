import numpy as np

# --- DADOS EXTRAÍDOS DOS CARTÕES DE VOO ---

# Voo 1 (CG Dianteiro ~20.15%) - Fonte: Pág 11 do Cartão
# Anotação "Pulso comandos LIVRES": tempos 0, 0.78, 2.2
# Interpretamos como: t0=0, t_pico1=0.78, t_pico2=2.2
t_voo1 = np.array([0.78, 2.2]) 
n_picos_voo1 = 2 # Piloto parou de anotar, oscilação morreu rápido

# Voo 2 (CG Traseiro ~22.99%) - Fonte: Pág 13 do Cartão
# Anotação "Doublet comandos LIVRES": tempos 0, 2.21, 3.58
# Interpretamos como: t_pico1=2.21, t_pico2=3.58
t_voo2 = np.array([2.21, 3.58])
n_picos_voo2 = 2 # Idem

# --- CÁLCULOS ---

def analisar_dinamica(tempos, n_picos):
    # 1. Período Amortecido (Td)
    # Diferença entre picos consecutivos
    Td = tempos[1] - tempos[0]
    
    # 2. Frequência Amortecida (wd)
    wd = (2 * np.pi) / Td
    
    # 3. Amortecimento (Zeta) - Método Contagem de Picos (PDF Eq. 36)
    # Se n=2, zeta = (7-2)/10 = 0.5. Se n=3, zeta=0.4.
    # Curto período é "heavily damped", então n é baixo.
    zeta = (7 - n_picos) / 10.0
    
    # Limite de segurança do método (zeta não pode ser > 1 nem < 0)
    zeta = max(0.1, min(zeta, 0.9))
    
    # 4. Frequência Natural (wn)
    wn = wd / np.sqrt(1 - zeta**2)
    
    return Td, wn, zeta

Td1, wn1, zeta1 = analisar_dinamica(t_voo1, n_picos_voo1)
Td2, wn2, zeta2 = analisar_dinamica(t_voo2, n_picos_voo2)

# --- OUTPUT ---
print("-" * 40)
print("ANÁLISE EXPERIMENTAL - CURTO PERÍODO")
print("-" * 40)
print(f"VOO 1 (CG Dianteiro):")
print(f"  Tempos entre picos: {Td1:.2f} s")
print(f"  Frequência Natural (wn): {wn1:.2f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta1:.2f}")
print("-" * 40)
print(f"VOO 2 (CG Traseiro):")
print(f"  Tempos entre picos: {Td2:.2f} s")
print(f"  Frequência Natural (wn): {wn2:.2f} rad/s")
print(f"  Amortecimento Est. (zeta): {zeta2:.2f}")
print("-" * 40)
print(f"Variação de Frequência: {((wn2-wn1)/wn1)*100:.1f}%")